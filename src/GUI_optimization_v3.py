"""
HPGe Monte Carlo Simulation GUI - Enhanced Multi-Source Version
Support for:
- Multi-isotopes mixing
- Custom multi-energy sources
- Automatic config.txt generation
- Professional ISOCS-style analysis
- Mode A/B Efficiency-based COI Analysis (Refactored)
- Apex-style live spectrum display during batch processing
"""

import sys
import json
import subprocess
import os
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
try:
    from scipy.integrate import trapezoid as trapz
except ImportError:
    from scipy.integrate import trapz
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                            QHBoxLayout, QPushButton, QLabel, QLineEdit, 
                            QComboBox, QGroupBox, QTabWidget, QTextEdit,
                            QTableWidget, QTableWidgetItem, QSpinBox, 
                            QDoubleSpinBox, QFileDialog, QProgressBar,
                            QMessageBox, QHeaderView, QDialog, QCheckBox,
                            QListWidget, QListWidgetItem, QGridLayout, 
                            QButtonGroup, QRadioButton, QSplitter, QScrollArea,
                            QProgressDialog, QInputDialog, QDialogButtonBox, QFrame)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont, QColor

# ==================================================================
# HELPER CLASSES
# ==================================================================

class ScrollableWidget(QWidget):
    """Wrapper to make any widget scrollable"""
    def __init__(self, content_widget):
        super().__init__()
        self.init_ui(content_widget)
    
    def init_ui(self, content_widget):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setWidget(content_widget)
        layout.addWidget(scroll)

class EfficiencyPlotDialog(QDialog):
    """Dialog to plot efficiency curve"""
    def __init__(self, peak_data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Efficiency vs Energy")
        self.resize(800, 600)
        self.peak_data = peak_data
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.plot_data()
        self.setLayout(layout)
        
    def plot_data(self):
        if not self.peak_data: return
        energies = [d['energy'] for d in self.peak_data]
        efficiencies = [d['efficiency'] for d in self.peak_data]
        names = [d['name'] for d in self.peak_data]
        
        self.ax.scatter(energies, efficiencies, color='red', s=50, zorder=5)
        for e, eff, name in zip(energies, efficiencies, names):
            self.ax.annotate(name, (e, eff), xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        if len(energies) >= 3:
            try:
                sorted_pairs = sorted(zip(energies, efficiencies))
                x_fit = np.array([p[0] for p in sorted_pairs])
                y_fit = np.array([p[1] for p in sorted_pairs])
                log_x = np.log(x_fit)
                log_y = np.log(y_fit)
                coeffs = np.polyfit(log_x, log_y, 2)
                poly = np.poly1d(coeffs)
                x_smooth = np.linspace(min(x_fit), max(x_fit), 100)
                y_smooth = np.exp(poly(np.log(x_smooth)))
                self.ax.plot(x_smooth, y_smooth, 'b--', alpha=0.7, label='Fit')
            except: pass

        self.ax.set_xlabel('Energy (keV)')
        self.ax.set_ylabel('Efficiency (Counts/Emission)')
        self.ax.set_title('Absolute Full-Energy Peak Efficiency')
        self.ax.grid(True, which='both', linestyle='--', alpha=0.5)
        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.canvas.draw()


# ==================================================================
# EMBEDDED LIVE SPECTRUM WIDGET (Apex-style)
# ==================================================================

class EmbeddedLiveSpectrum(QWidget):
    """Compact embedded live spectrum display for batch processing"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.energies = []
        self.total_counts = 0
        self.start_time = None
        self.is_running = False
        self.output_file = None
        self.last_file_position = 0
        self.update_interval = 500
        
        # Spectrum parameters
        self.num_bins = 1024
        self.energy_max = 2500
        self.bin_width = self.energy_max / self.num_bins
        self.spectrum = np.zeros(self.num_bins)
        
        # Current isotope/mode info
        self.current_isotope = ""
        self.current_mode = ""
        self.batch_progress = (0, 0)  # (current, total)
        
        self.init_ui()
        
        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_spectrum)
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(5)
        
        # Status bar with isotope and mode info
        status_frame = QFrame()
        status_frame.setStyleSheet("""
            QFrame {
                background-color: #1a1a2e;
                border-radius: 5px;
                padding: 5px;
            }
        """)
        status_layout = QHBoxLayout(status_frame)
        status_layout.setContentsMargins(10, 5, 10, 5)
        
        # Status indicator
        self.status_indicator = QLabel("⏹")
        self.status_indicator.setStyleSheet("font-size: 16px; color: #666;")
        status_layout.addWidget(self.status_indicator)
        
        # Current isotope label
        self.isotope_label = QLabel("No simulation running")
        self.isotope_label.setStyleSheet("font-weight: bold; font-size: 14px; color: #00ff00;")
        status_layout.addWidget(self.isotope_label)
        
        status_layout.addStretch()
        
        # Mode indicator (A or B)
        self.mode_label = QLabel("")
        self.mode_label.setStyleSheet("font-weight: bold; font-size: 14px; color: #ffaa00; padding: 2px 8px; background-color: #333; border-radius: 3px;")
        status_layout.addWidget(self.mode_label)
        
        # Progress indicator
        self.progress_label = QLabel("")
        self.progress_label.setStyleSheet("font-size: 12px; color: #aaa;")
        status_layout.addWidget(self.progress_label)
        
        layout.addWidget(status_frame)
        
        # Stats bar
        stats_frame = QFrame()
        stats_frame.setStyleSheet("background-color: #16213e; border-radius: 3px; padding: 3px;")
        stats_layout = QHBoxLayout(stats_frame)
        stats_layout.setContentsMargins(10, 3, 10, 3)
        
        self.stats_label = QLabel("Counts: 0 | Rate: 0 cps | Time: 0s")
        self.stats_label.setStyleSheet("font-family: monospace; font-size: 11px; color: #00ff00;")
        stats_layout.addWidget(self.stats_label)
        
        stats_layout.addStretch()
        
        # Y-scale toggle
        self.log_scale_cb = QCheckBox("Log Y")
        self.log_scale_cb.setStyleSheet("color: #aaa; font-size: 10px;")
        self.log_scale_cb.stateChanged.connect(self.update_plot_only)
        stats_layout.addWidget(self.log_scale_cb)
        
        layout.addWidget(stats_frame)
        
        # Matplotlib figure with compact dark theme
        self.figure, self.ax = plt.subplots(figsize=(10, 3))
        self.figure.patch.set_facecolor('#1a1a2e')
        self.ax.set_facecolor('#16213e')
        self.figure.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.15)
        
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumHeight(200)
        self.canvas.setMaximumHeight(250)
        layout.addWidget(self.canvas)
        
        self.init_plot()
        self.setVisible(False)  # Hidden by default
    
    def init_plot(self):
        """Initialize the spectrum plot"""
        self.ax.clear()
        
        self.energy_bins = np.linspace(0, self.energy_max, self.num_bins + 1)
        self.energy_centers = (self.energy_bins[:-1] + self.energy_bins[1:]) / 2
        
        self.line, = self.ax.step(self.energy_centers, self.spectrum, 
                                   where='mid', color='#00ff00', linewidth=0.8)
        
        self.fill = self.ax.fill_between(self.energy_centers, 0, self.spectrum,
                                          step='mid', alpha=0.3, color='#00ff00')
        
        self.ax.set_xlabel('Energy (keV)', color='white', fontsize=9)
        self.ax.set_ylabel('Counts', color='white', fontsize=9)
        self.ax.tick_params(colors='white', labelsize=8)
        self.ax.spines['bottom'].set_color('white')
        self.ax.spines['top'].set_color('white')
        self.ax.spines['left'].set_color('white')
        self.ax.spines['right'].set_color('white')
        self.ax.grid(True, alpha=0.3, color='gray', linestyle='--')
        
        self.ax.set_xlim(0, 1500)
        self.ax.set_ylim(0, 100)
        
        self.canvas.draw()
    
    def set_current_simulation(self, isotope, mode, current_idx, total):
        """Update the display with current simulation info"""
        self.current_isotope = isotope
        self.current_mode = mode
        self.batch_progress = (current_idx, total)
        
        self.isotope_label.setText(f"📊 {isotope}")
        self.mode_label.setText(f"Mode {mode}")
        
        if mode == "A":
            self.mode_label.setStyleSheet("font-weight: bold; font-size: 14px; color: #00ff00; padding: 2px 8px; background-color: #1a4d1a; border-radius: 3px;")
        else:
            self.mode_label.setStyleSheet("font-weight: bold; font-size: 14px; color: #00aaff; padding: 2px 8px; background-color: #1a3d5c; border-radius: 3px;")
        
        self.progress_label.setText(f"[{current_idx}/{total}]")
    
    def start_live_update(self, output_file):
        """Start live spectrum updates"""
        self.output_file = output_file
        self.is_running = True
        self.start_time = time.time()
        self.last_file_position = 0
        self.spectrum = np.zeros(self.num_bins)
        self.total_counts = 0
        self.energies = []
        
        self.status_indicator.setText("🔴")
        self.status_indicator.setStyleSheet("font-size: 16px; color: #ff0000;")
        
        self.setVisible(True)
        self.update_timer.start(self.update_interval)
    
    def stop_live_update(self):
        """Stop live updates"""
        self.is_running = False
        self.update_timer.stop()
        
        self.status_indicator.setText("✓")
        self.status_indicator.setStyleSheet("font-size: 16px; color: #00ff00;")
        
        self.update_spectrum()
    
    def finish_batch(self):
        """Called when batch is complete"""
        self.isotope_label.setText("✅ Batch Complete!")
        self.mode_label.setText("")
        self.status_indicator.setText("✓")
        self.status_indicator.setStyleSheet("font-size: 16px; color: #00ff00;")
    
    def update_spectrum(self):
        """Read new data and update spectrum"""
        if not self.output_file or not os.path.exists(self.output_file):
            return
        
        try:
            with open(self.output_file, 'r') as f:
                f.seek(self.last_file_position)
                new_lines = f.readlines()
                self.last_file_position = f.tell()
            
            new_energies = []
            for line in new_lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                try:
                    parts = line.split()
                    if len(parts) >= 3:
                        n_det = int(float(parts[2]))
                        for i in range(n_det):
                            idx = 4 + i * 3
                            if idx < len(parts):
                                energy = float(parts[idx])
                                new_energies.append(energy)
                except:
                    continue
            
            if new_energies:
                for e in new_energies:
                    bin_idx = int(e / self.bin_width)
                    if 0 <= bin_idx < self.num_bins:
                        self.spectrum[bin_idx] += 1
                        self.total_counts += 1
                
                self.energies.extend(new_energies)
            
            self.update_plot_only()
            
            elapsed = time.time() - self.start_time if self.start_time else 0
            rate = self.total_counts / elapsed if elapsed > 0 else 0
            self.stats_label.setText(f"Counts: {self.total_counts:,} | Rate: {rate:.1f} cps | Time: {elapsed:.1f}s")
            
        except Exception as e:
            print(f"Update error: {e}")
    
    def update_plot_only(self):
        """Update plot without reading file"""
        self.line.set_ydata(self.spectrum)
        
        self.fill.remove()
        self.fill = self.ax.fill_between(self.energy_centers, 0, self.spectrum,
                                          step='mid', alpha=0.3, color='#00ff00')
        
        if self.log_scale_cb.isChecked():
            self.ax.set_yscale('log')
            self.ax.set_ylim(0.1, max(10, self.spectrum.max() * 1.5))
        else:
            self.ax.set_yscale('linear')
            visible_max = self.spectrum.max()
            self.ax.set_ylim(0, max(10, visible_max * 1.1))
        
        self.canvas.draw_idle()
    
    def clear_spectrum(self):
        """Clear the spectrum for new isotope"""
        self.spectrum = np.zeros(self.num_bins)
        self.total_counts = 0
        self.energies = []
        self.last_file_position = 0
        self.update_plot_only()
        self.stats_label.setText("Counts: 0 | Rate: 0 cps | Time: 0s")


# ==================================================================
# BATCH RESULTS DIALOG (Enhanced)
# ==================================================================

class BatchResultsDialog(QDialog):
    """Ultra-Professional Batch Results Dialog with Plotting and Export"""
    
    def __init__(self, parent, results_data, beamon):
        super().__init__(parent)
        self.setWindowTitle("📊 Batch COI Analysis Results")
        
        self.setWindowFlags(Qt.Window | Qt.WindowMinimizeButtonHint | 
                           Qt.WindowMaximizeButtonHint | Qt.WindowCloseButtonHint |
                           Qt.WindowSystemMenuHint)
        
        self.setSizeGripEnabled(True)
        self.setMinimumSize(800, 600)
        self.resize(1200, 800)
        
        self.results_data = results_data
        self.beamon = beamon
        self.experimental_data = []
        self.cascade_isotopes = ['Co60', 'Y88', 'Na22', 'Eu152', 'Ba133']
        
        self.processed_data = self.process_results()
        self.auto_save_results()
        self.init_ui()
    
    def auto_save_results(self):
        """Auto-save results to JSON file for later retrieval"""
        try:
            save_data = {
                'beamon': self.beamon,
                'results': self.results_data,
                'processed': self.processed_data,
                'timestamp': str(np.datetime64('now'))
            }
            
            save_path = os.path.join(os.getcwd(), 'batch_results_cache.json')
            
            def convert_to_serializable(obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                elif isinstance(obj, dict):
                    return {k: convert_to_serializable(v) for k, v in obj.items()}
                elif isinstance(obj, list):
                    return [convert_to_serializable(i) for i in obj]
                return obj
            
            save_data = convert_to_serializable(save_data)
            
            with open(save_path, 'w') as f:
                json.dump(save_data, f, indent=2)
            
            print(f"Results auto-saved to: {save_path}")
        except Exception as e:
            print(f"Warning: Could not auto-save results: {e}")
    
    def process_results(self):
        """Process raw results into structured data with calculations"""
        processed = []
        
        for isotope, data in self.results_data.items():
            peaks = ISOTOPE_DATABASE.get(isotope, [])
            num_energies = len(peaks)
            has_cascade = isotope in self.cascade_isotopes
            
            for energy_data in data['energies']:
                energy = energy_data['energy']
                intensity = energy_data['intensity']
                cnt_a = energy_data['cnt_a']
                cnt_b = energy_data['cnt_b']
                
                if num_energies == 1:
                    emissions = self.beamon
                else:
                    emissions = self.beamon * (intensity / 100.0)
                
                # Mode B correction: Custom mode divides beamon among energies
                # So cnt_b must be multiplied by num_energies for correct efficiency
                if num_energies > 1:
                    cnt_b_adj = cnt_b * num_energies
                else:
                    cnt_b_adj = cnt_b
                
                eff_a = cnt_a / emissions if emissions > 0 else 0
                eff_b = cnt_b_adj / emissions if emissions > 0 else 0  # Use adjusted count!
                sigma_eff_a = eff_a * np.sqrt(1.0 / cnt_a) if cnt_a > 0 else 0
                sigma_eff_b = eff_b * np.sqrt(1.0 / cnt_b_adj) if cnt_b_adj > 0 else 0
                
                coi = cnt_a / cnt_b_adj if cnt_b_adj > 0 else 1.0
                
                if cnt_a > 0 and cnt_b_adj > 0:
                    sigma_coi = coi * np.sqrt(1.0/cnt_a + 1.0/cnt_b_adj)
                else:
                    sigma_coi = 0
                
                processed.append({
                    'isotope': isotope,
                    'energy': energy,
                    'intensity': intensity,
                    'cnt_a': cnt_a,
                    'cnt_b': cnt_b,
                    'emissions': emissions,
                    'eff_a': eff_a,
                    'eff_b': eff_b,
                    'sigma_eff_a': sigma_eff_a,
                    'sigma_eff_b': sigma_eff_b,
                    'coi': coi,
                    'sigma_coi': sigma_coi,
                    'has_cascade': has_cascade,
                    'num_energies': num_energies
                })
        
        processed.sort(key=lambda x: x['energy'])
        return processed
    
    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(5, 5, 5, 5)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        
        content_widget = QWidget()
        layout = QVBoxLayout(content_widget)
        
        title = QLabel("<h2>📊 Multi-Isotope Batch Analysis Results</h2>")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        summary = QLabel(f"<b>Isotopes:</b> {len(self.results_data)} | "
                        f"<b>Energy Lines:</b> {len(self.processed_data)} | "
                        f"<b>BeamOn:</b> {self.beamon:,}")
        summary.setAlignment(Qt.AlignCenter)
        summary.setStyleSheet("background-color: #E3F2FD; padding: 10px; border-radius: 5px;")
        layout.addWidget(summary)
        
        self.tabs = QTabWidget()
        self.tabs.setMinimumHeight(500)
        
        self.tabs.addTab(self.create_table_tab(), "📋 Results Table")
        self.tabs.addTab(self.create_efficiency_plot_tab(), "📈 Efficiency Plot")
        self.tabs.addTab(self.create_coi_plot_tab(), "📉 COI Plot")
        self.tabs.addTab(self.create_comparison_tab(), "🔬 Exp. Comparison")
        
        layout.addWidget(self.tabs)
        
        scroll_area.setWidget(content_widget)
        main_layout.addWidget(scroll_area)
        
        btn_layout = QHBoxLayout()
        
        export_csv_btn = QPushButton("📄 Export CSV")
        export_csv_btn.clicked.connect(self.export_csv)
        #export_csv_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; padding: 8px;")
        btn_layout.addWidget(export_csv_btn)
        
        save_json_btn = QPushButton("💾 Save Results (JSON)")
        save_json_btn.clicked.connect(self.save_results_json)
        #save_json_btn.setStyleSheet("background-color: #9C27B0; color: white; font-weight: bold; padding: 8px;")
        btn_layout.addWidget(save_json_btn)
        
        export_plot_btn = QPushButton("🖼️ Export Plot")
        export_plot_btn.clicked.connect(self.export_plot)
        #export_plot_btn.setStyleSheet("background-color: #2196F3; color: white; font-weight: bold; padding: 8px;")
        btn_layout.addWidget(export_plot_btn)
        
        import_exp_btn = QPushButton("📥 Import Experimental Data")
        import_exp_btn.clicked.connect(self.import_experimental)
        #import_exp_btn.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold; padding: 8px;")
        btn_layout.addWidget(import_exp_btn)
        
        close_btn = QPushButton("✖ Close")
        close_btn.clicked.connect(self.accept)
        #close_btn.setStyleSheet("background-color: #757575; color: white; font-weight: bold; padding: 8px;")
        btn_layout.addWidget(close_btn)
        
        main_layout.addLayout(btn_layout)
        self.setLayout(main_layout)
    
    def create_table_tab(self):
        """Create the results table tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        self.table = QTableWidget()
        self.table.setColumnCount(12)
        self.table.setHorizontalHeaderLabels([
            "Isotope", "Energy\n(keV)", "Intensity\n(%)", "Emissions",
            "Counts A", "Counts B", 
            "Eff A", "σ(Eff A)", "Eff B", "σ(Eff B)",
            "COI", "COI ± σ"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.setAlternatingRowColors(True)
        self.table.setStyleSheet("QTableWidget { gridline-color: #d0d0d0; }"
                                  "QTableWidget::item:alternate { background-color: #f5f5f5; }")
        
        for row_idx, d in enumerate(self.processed_data):
            self.table.insertRow(row_idx)
            
            self.table.setItem(row_idx, 0, QTableWidgetItem(d['isotope']))
            self.table.setItem(row_idx, 1, QTableWidgetItem(f"{d['energy']:.3f}"))
            self.table.setItem(row_idx, 2, QTableWidgetItem(f"{d['intensity']:.2f}"))
            self.table.setItem(row_idx, 3, QTableWidgetItem(f"{int(d['emissions']):,}"))
            self.table.setItem(row_idx, 4, QTableWidgetItem(f"{int(d['cnt_a']):,}"))
            self.table.setItem(row_idx, 5, QTableWidgetItem(f"{int(d['cnt_b']):,}"))
            self.table.setItem(row_idx, 6, QTableWidgetItem(f"{d['eff_a']:.6f}"))
            self.table.setItem(row_idx, 7, QTableWidgetItem(f"{d['sigma_eff_a']:.6f}"))
            self.table.setItem(row_idx, 8, QTableWidgetItem(f"{d['eff_b']:.6f}"))
            self.table.setItem(row_idx, 9, QTableWidgetItem(f"{d['sigma_eff_b']:.6f}"))
            
            coi_item = QTableWidgetItem(f"{d['coi']:.4f}")
            if d['has_cascade'] and abs(d['coi'] - 1.0) > 0.01:
                coi_item.setBackground(QColor("#fff3cd"))
            elif not d['has_cascade'] and abs(d['coi'] - 1.0) > 0.05:
                coi_item.setBackground(QColor("#f8d7da"))
            else:
                coi_item.setBackground(QColor("#d4edda"))
            self.table.setItem(row_idx, 10, coi_item)
            
            self.table.setItem(row_idx, 11, QTableWidgetItem(f"{d['coi']:.4f} ± {d['sigma_coi']:.4f}"))
        
        layout.addWidget(self.table)
        widget.setLayout(layout)
        return widget
    
    def create_efficiency_plot_tab(self):
        """Create the efficiency plot tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        ctrl_layout = QHBoxLayout()
        
        ctrl_layout.addWidget(QLabel("Plot Type:"))
        self.eff_plot_type = QComboBox()
        self.eff_plot_type.addItems(["Mode A Only", "Mode B Only", "Both Modes", "Mode A with Fit"])
        self.eff_plot_type.currentIndexChanged.connect(self.update_efficiency_plot)
        ctrl_layout.addWidget(self.eff_plot_type)
        
        ctrl_layout.addWidget(QLabel("Scale:"))
        self.eff_scale = QComboBox()
        self.eff_scale.addItems(["Log-Log", "Linear", "Log-Linear", "Linear-Log"])
        self.eff_scale.currentIndexChanged.connect(self.update_efficiency_plot)
        ctrl_layout.addWidget(self.eff_scale)
        
        self.show_labels_cb = QCheckBox("Show Labels")
        self.show_labels_cb.setChecked(True)
        self.show_labels_cb.stateChanged.connect(self.update_efficiency_plot)
        ctrl_layout.addWidget(self.show_labels_cb)
        
        self.show_error_bars_cb = QCheckBox("Show Error Bars")
        self.show_error_bars_cb.setChecked(True)
        self.show_error_bars_cb.stateChanged.connect(self.update_efficiency_plot)
        ctrl_layout.addWidget(self.show_error_bars_cb)
        
        ctrl_layout.addStretch()
        layout.addLayout(ctrl_layout)
        
        self.eff_figure, self.eff_ax = plt.subplots(figsize=(10, 6))
        self.eff_canvas = FigureCanvasQTAgg(self.eff_figure)
        self.eff_toolbar = NavigationToolbar2QT(self.eff_canvas, widget)
        layout.addWidget(self.eff_toolbar)
        layout.addWidget(self.eff_canvas)
        
        self.update_efficiency_plot()
        
        widget.setLayout(layout)
        return widget
    
    def update_efficiency_plot(self):
        """Update the efficiency plot"""
        self.eff_ax.clear()
        
        energies = [d['energy'] for d in self.processed_data]
        eff_a = [d['eff_a'] for d in self.processed_data]
        eff_b = [d['eff_b'] for d in self.processed_data]
        sigma_a = [d['sigma_eff_a'] for d in self.processed_data]
        sigma_b = [d['sigma_eff_b'] for d in self.processed_data]
        labels = [d['isotope'] for d in self.processed_data]
        
        plot_type = self.eff_plot_type.currentText()
        show_errors = self.show_error_bars_cb.isChecked()
        
        if plot_type in ["Mode A Only", "Both Modes", "Mode A with Fit"]:
            if show_errors:
                self.eff_ax.errorbar(energies, eff_a, yerr=sigma_a, fmt='o', 
                                     color='#1976D2', markersize=8, capsize=3,
                                     label='Mode A (Simulation)', zorder=5)
            else:
                self.eff_ax.scatter(energies, eff_a, color='#1976D2', s=60, 
                                   label='Mode A (Simulation)', zorder=5)
        
        if plot_type in ["Mode B Only", "Both Modes"]:
            if show_errors:
                self.eff_ax.errorbar(energies, eff_b, yerr=sigma_b, fmt='s', 
                                     color='#388E3C', markersize=8, capsize=3,
                                     label='Mode B (Independent)', zorder=5)
            else:
                self.eff_ax.scatter(energies, eff_b, color='#388E3C', s=60, marker='s',
                                   label='Mode B (Independent)', zorder=5)
        
        if self.experimental_data:
            exp_e = [d['energy'] for d in self.experimental_data]
            exp_eff = [d['efficiency'] for d in self.experimental_data]
            exp_sigma = [d.get('sigma', 0) for d in self.experimental_data]
            
            if show_errors and any(exp_sigma):
                self.eff_ax.errorbar(exp_e, exp_eff, yerr=exp_sigma, fmt='^', 
                                     color='#D32F2F', markersize=10, capsize=3,
                                     label='Experimental', zorder=6)
            else:
                self.eff_ax.scatter(exp_e, exp_eff, color='#D32F2F', s=80, marker='^',
                                   label='Experimental', zorder=6)
        
        if plot_type == "Mode A with Fit" and len(energies) >= 3:
            try:
                sorted_data = sorted(zip(energies, eff_a))
                x_fit = np.array([p[0] for p in sorted_data])
                y_fit = np.array([p[1] for p in sorted_data])
                
                log_x = np.log(x_fit)
                log_y = np.log(y_fit)
                coeffs = np.polyfit(log_x, log_y, 3)
                poly = np.poly1d(coeffs)
                
                x_smooth = np.logspace(np.log10(min(x_fit)*0.9), np.log10(max(x_fit)*1.1), 200)
                y_smooth = np.exp(poly(np.log(x_smooth)))
                
                self.eff_ax.plot(x_smooth, y_smooth, '--', color='#1976D2', alpha=0.7, 
                                linewidth=2, label='Polynomial Fit')
            except Exception as e:
                print(f"Fit error: {e}")
        
        if self.show_labels_cb.isChecked():
            for e, eff, lbl in zip(energies, eff_a if "Mode A" in plot_type or plot_type == "Both Modes" else eff_b, labels):
                self.eff_ax.annotate(lbl, (e, eff), xytext=(5, 5), 
                                    textcoords='offset points', fontsize=8, alpha=0.8)
        
        scale = self.eff_scale.currentText()
        if scale == "Log-Log":
            self.eff_ax.set_xscale('log')
            self.eff_ax.set_yscale('log')
        elif scale == "Log-Linear":
            self.eff_ax.set_xscale('log')
            self.eff_ax.set_yscale('linear')
        elif scale == "Linear-Log":
            self.eff_ax.set_xscale('linear')
            self.eff_ax.set_yscale('log')
        else:
            self.eff_ax.set_xscale('linear')
            self.eff_ax.set_yscale('linear')
        
        self.eff_ax.set_xlabel('Energy (keV)', fontsize=12)
        self.eff_ax.set_ylabel('Efficiency', fontsize=12)
        self.eff_ax.set_title('Full-Energy Peak Efficiency', fontsize=14, fontweight='bold')
        self.eff_ax.grid(True, which='both', linestyle='--', alpha=0.5)
        self.eff_ax.legend(loc='best')
        
        self.eff_figure.tight_layout()
        self.eff_canvas.draw()
    
    def create_coi_plot_tab(self):
        """Create the COI plot tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        self.coi_figure, self.coi_ax = plt.subplots(figsize=(10, 6))
        self.coi_canvas = FigureCanvasQTAgg(self.coi_figure)
        self.coi_toolbar = NavigationToolbar2QT(self.coi_canvas, widget)
        layout.addWidget(self.coi_toolbar)
        layout.addWidget(self.coi_canvas)
        
        self.update_coi_plot()
        
        widget.setLayout(layout)
        return widget
    
    def update_coi_plot(self):
        """Update the COI plot"""
        self.coi_ax.clear()
        
        energies = [d['energy'] for d in self.processed_data]
        coi_values = [d['coi'] for d in self.processed_data]
        sigma_coi = [d['sigma_coi'] for d in self.processed_data]
        isotopes = [d['isotope'] for d in self.processed_data]
        has_cascade = [d['has_cascade'] for d in self.processed_data]
        
        colors = ['#FF5722' if hc else '#2196F3' for hc in has_cascade]
        
        self.coi_ax.errorbar(energies, coi_values, yerr=sigma_coi, fmt='o', 
                            color='#333', markersize=0, capsize=3, zorder=4)
        
        scatter = self.coi_ax.scatter(energies, coi_values, c=colors, s=80, zorder=5, edgecolors='black')
        
        self.coi_ax.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label='COI = 1.0 (No cascade)')
        
        for e, coi, lbl in zip(energies, coi_values, isotopes):
            self.coi_ax.annotate(lbl, (e, coi), xytext=(5, 5), 
                                textcoords='offset points', fontsize=8, alpha=0.8)
        
        self.coi_ax.set_xlabel('Energy (keV)', fontsize=12)
        self.coi_ax.set_ylabel('COI Factor', fontsize=12)
        self.coi_ax.set_title('Coincidence Summing Correction Factor', fontsize=14, fontweight='bold')
        self.coi_ax.grid(True, which='both', linestyle='--', alpha=0.5)
        
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#FF5722', edgecolor='black', label='Cascade isotopes'),
            Patch(facecolor='#2196F3', edgecolor='black', label='Single-gamma isotopes'),
        ]
        self.coi_ax.legend(handles=legend_elements, loc='best')
        
        self.coi_figure.tight_layout()
        self.coi_canvas.draw()
    
    def create_comparison_tab(self):
        """Create the experimental comparison tab"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        info_label = QLabel("Import experimental data to compare with simulation results.")
        info_label.setStyleSheet("font-style: italic; color: #666;")
        layout.addWidget(info_label)
        
        self.exp_table = QTableWidget()
        self.exp_table.setColumnCount(4)
        self.exp_table.setHorizontalHeaderLabels(["Energy (keV)", "Exp. Efficiency", "σ", "Sim. Efficiency"])
        self.exp_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        layout.addWidget(self.exp_table)
        
        self.comp_figure, self.comp_ax = plt.subplots(figsize=(10, 5))
        self.comp_canvas = FigureCanvasQTAgg(self.comp_figure)
        layout.addWidget(self.comp_canvas)
        
        self.stats_label = QLabel("")
        self.stats_label.setStyleSheet("padding: 10px; background-color: #f5f5f5; border-radius: 5px;")
        layout.addWidget(self.stats_label)
        
        widget.setLayout(layout)
        return widget
    
    def import_experimental(self):
        """Import experimental efficiency data"""
        filename, _ = QFileDialog.getOpenFileName(self, "Import Experimental Data", 
                                                   "", "CSV Files (*.csv);;Text Files (*.txt);;All Files (*)")
        if not filename:
            return
        
        try:
            self.experimental_data = []
            
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.replace(',', ' ').split()
                if len(parts) >= 2:
                    try:
                        energy = float(parts[0])
                        efficiency = float(parts[1])
                        sigma = float(parts[2]) if len(parts) >= 3 else 0
                        self.experimental_data.append({
                            'energy': energy,
                            'efficiency': efficiency,
                            'sigma': sigma
                        })
                    except ValueError:
                        continue
            
            if self.experimental_data:
                QMessageBox.information(self, "Import Success", 
                                       f"Imported {len(self.experimental_data)} data points.")
                self.update_efficiency_plot()
                self.update_comparison_tab()
            else:
                QMessageBox.warning(self, "Import Failed", "No valid data found in file.")
                
        except Exception as e:
            QMessageBox.critical(self, "Import Error", f"Error reading file:\n{str(e)}")
    
    def update_comparison_tab(self):
        """Update comparison tab with experimental data"""
        self.exp_table.setRowCount(0)
        
        sim_efficiencies = {d['energy']: d['eff_a'] for d in self.processed_data}
        residuals = []
        
        for row_idx, exp in enumerate(self.experimental_data):
            self.exp_table.insertRow(row_idx)
            
            energy = exp['energy']
            exp_eff = exp['efficiency']
            exp_sigma = exp.get('sigma', 0)
            
            closest_sim_e = min(sim_efficiencies.keys(), key=lambda x: abs(x - energy))
            sim_eff = sim_efficiencies[closest_sim_e] if abs(closest_sim_e - energy) < 5 else None
            
            self.exp_table.setItem(row_idx, 0, QTableWidgetItem(f"{energy:.3f}"))
            self.exp_table.setItem(row_idx, 1, QTableWidgetItem(f"{exp_eff:.6f}"))
            self.exp_table.setItem(row_idx, 2, QTableWidgetItem(f"{exp_sigma:.6f}" if exp_sigma else "N/A"))
            
            if sim_eff is not None:
                self.exp_table.setItem(row_idx, 3, QTableWidgetItem(f"{sim_eff:.6f}"))
                residuals.append((exp_eff - sim_eff) / exp_eff * 100 if exp_eff != 0 else 0)
            else:
                self.exp_table.setItem(row_idx, 3, QTableWidgetItem("No match"))
        
        self.comp_ax.clear()
        
        exp_e = [d['energy'] for d in self.experimental_data]
        exp_eff = [d['efficiency'] for d in self.experimental_data]
        exp_sigma = [d.get('sigma', 0) for d in self.experimental_data]
        
        sim_e = [d['energy'] for d in self.processed_data]
        sim_eff = [d['eff_a'] for d in self.processed_data]
        sim_sigma = [d['sigma_eff_a'] for d in self.processed_data]
        
        self.comp_ax.errorbar(exp_e, exp_eff, yerr=exp_sigma if any(exp_sigma) else None, 
                             fmt='^', color='#D32F2F', markersize=10, capsize=4,
                             label='Experimental', zorder=6)
        
        self.comp_ax.errorbar(sim_e, sim_eff, yerr=sim_sigma, 
                             fmt='o', color='#1976D2', markersize=8, capsize=3,
                             label='Simulation (Mode A)', zorder=5)
        
        self.comp_ax.set_xlabel('Energy (keV)', fontsize=12)
        self.comp_ax.set_ylabel('Efficiency', fontsize=12)
        self.comp_ax.set_title('Simulation vs Experimental Efficiency', fontsize=14, fontweight='bold')
        self.comp_ax.set_xscale('log')
        self.comp_ax.set_yscale('log')
        self.comp_ax.grid(True, which='both', linestyle='--', alpha=0.5)
        self.comp_ax.legend(loc='best')
        
        self.comp_figure.tight_layout()
        self.comp_canvas.draw()
        
        if residuals:
            mean_res = np.mean(residuals)
            std_res = np.std(residuals)
            self.stats_label.setText(
                f"<b>Comparison Statistics:</b><br>"
                f"Mean Residual: {mean_res:.2f}%<br>"
                f"Std Residual: {std_res:.2f}%<br>"
                f"Data Points Matched: {len(residuals)}/{len(self.experimental_data)}"
            )
    
    def export_csv(self):
        """Export results to CSV"""
        filename, _ = QFileDialog.getSaveFileName(self, "Export Results", 
                                                   "batch_results.csv", "CSV (*.csv)")
        if filename:
            with open(filename, 'w') as f:
                f.write("Isotope,Energy_keV,Intensity_pct,Emissions,Counts_A,Counts_B,"
                       "Eff_A,Sigma_Eff_A,Eff_B,Sigma_Eff_B,COI,Sigma_COI,Has_Cascade\n")
                
                for d in self.processed_data:
                    f.write(f"{d['isotope']},{d['energy']:.3f},{d['intensity']:.2f},"
                           f"{int(d['emissions'])},{int(d['cnt_a'])},{int(d['cnt_b'])},"
                           f"{d['eff_a']:.8f},{d['sigma_eff_a']:.8f},"
                           f"{d['eff_b']:.8f},{d['sigma_eff_b']:.8f},"
                           f"{d['coi']:.6f},{d['sigma_coi']:.6f},{d['has_cascade']}\n")
            
            QMessageBox.information(self, "Export", f"Results exported to:\n{filename}")

    def save_results_json(self):
        """Save results to JSON file for later loading"""
        filename, _ = QFileDialog.getSaveFileName(self, "Save Results", 
                                                   "batch_results.json", "JSON (*.json)")
        if filename:
            try:
                save_data = {
                    'beamon': self.beamon,
                    'results': self.results_data,
                    'processed': self.processed_data,
                    'timestamp': str(np.datetime64('now'))
                }
                
                def convert_to_serializable(obj):
                    if isinstance(obj, np.integer):
                        return int(obj)
                    elif isinstance(obj, np.floating):
                        return float(obj)
                    elif isinstance(obj, np.ndarray):
                        return obj.tolist()
                    elif isinstance(obj, dict):
                        return {k: convert_to_serializable(v) for k, v in obj.items()}
                    elif isinstance(obj, list):
                        return [convert_to_serializable(i) for i in obj]
                    return obj
                
                save_data = convert_to_serializable(save_data)
                
                with open(filename, 'w') as f:
                    json.dump(save_data, f, indent=2)
                
                QMessageBox.information(self, "Save", f"Results saved to:\n{filename}")
            except Exception as e:
                QMessageBox.critical(self, "Save Error", f"Error saving file:\n{str(e)}")
    
    def export_plot(self):
        """Export current plot to image"""
        current_tab = self.tabs.currentIndex()
        
        if current_tab == 1:
            figure = self.eff_figure
            default_name = "efficiency_plot.png"
        elif current_tab == 2:
            figure = self.coi_figure
            default_name = "coi_plot.png"
        elif current_tab == 3:
            figure = self.comp_figure
            default_name = "comparison_plot.png"
        else:
            QMessageBox.warning(self, "Export", "Please select a plot tab to export.")
            return
        
        filename, _ = QFileDialog.getSaveFileName(self, "Export Plot", 
                                                   default_name, 
                                                   "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if filename:
            figure.savefig(filename, dpi=300, bbox_inches='tight')
            QMessageBox.information(self, "Export", f"Plot exported to:\n{filename}")


# ==================================================================
# DATA & PHYSICS FUNCTIONS
# ==================================================================

ISOTOPE_DATABASE = {
    'Cs137': [{'energy': 661.657, 'intensity': 85.1, 'name': 'Cs-137'}],
    'Co60': [{'energy': 1173.228, 'intensity': 99.85, 'name': 'Co-60 (1)'},
             {'energy': 1332.492, 'intensity': 99.9826, 'name': 'Co-60 (2)'}],
    'Na22': [{'energy': 1274.537, 'intensity': 99.94, 'name': 'Na-22'}],
    'Am241': [{'energy': 59.5409, 'intensity': 35.9, 'name': 'Am-241'}],
    'Ba133': [{'energy': 30.625, 'intensity': 33.1, 'name': 'Ba-133 (Cs K)'},
              {'energy': 80.9979, 'intensity': 32.9, 'name': 'Ba-133'},
              {'energy': 356.017, 'intensity': 62.05, 'name': 'Ba-133'}],
    'Eu152': [{'energy': 121.782, 'intensity': 28.58, 'name': 'Eu-152'},
              {'energy': 344.278, 'intensity': 26.57, 'name': 'Eu-152'},
              {'energy': 1408.013, 'intensity': 21.01, 'name': 'Eu-152'}],
    'Mn54': [{'energy': 834.848, 'intensity': 99.976, 'name': 'Mn-54'}],
    'Cd109': [{'energy': 88.034, 'intensity': 3.644, 'name': 'Cd-109'}],
    'Sn113': [{'energy': 391.698, 'intensity': 64.97, 'name': 'Sn-113'}],
    'Ce139': [{'energy': 165.8575, 'intensity': 79.90, 'name': 'Ce-139'}],
    'Sr85': [{'energy': 514.0048, 'intensity': 96, 'name': 'Sr-85'}],
    'Zn65': [{'energy': 1115.52, 'intensity': 50.04, 'name': 'Zn-65'}],
    'Co57': [{'energy': 122.06065, 'intensity': 85.60, 'name': 'Co-57 (1)'},
             {'energy': 136.4736, 'intensity': 10.68, 'name': 'Co-57 (2)'}],
    'Y88': [{'energy': 898.042, 'intensity': 93.7, 'name': 'Y-88'},
            {'energy': 1836.063, 'intensity': 99.2, 'name': 'Y-88'}],
}

def calculate_peak_area(bins, counts, peak_energy, window_percent=0.03):
    window_low = peak_energy * (1 - window_percent)
    window_high = peak_energy * (1 + window_percent)
    mask = (bins >= window_low) & (bins <= window_high)
    if not np.any(mask): return 0, 0, 0
    
    window_bins = bins[mask]
    window_counts = counts[mask]
    if len(window_counts) < 3: return np.sum(window_counts), np.sum(window_counts), 0
    
    bg_left = np.mean(window_counts[:3])
    bg_right = np.mean(window_counts[-3:])
    background = np.linspace(bg_left, bg_right, len(window_counts))
    gross_area = trapz(window_counts, window_bins)
    bg_area = trapz(background, window_bins)
    return gross_area - bg_area, gross_area, bg_area


# ==================================================================
# PLOTTING CLASSES
# ==================================================================

class InteractiveSpectrum(FigureCanvasQTAgg):
    peak_selected = pyqtSignal(float, float)
    
    def __init__(self, parent=None):
        self.figure, self.ax = plt.subplots(figsize=(10, 6))
        super().__init__(self.figure)
        self.setParent(parent)
        self.energies, self.spectrum_data, self.bins = None, None, None
        self.detected_peaks = []
        self.selecting = False
        self.select_start = None
        self.select_rect = None
        
        self.mpl_connect('button_press_event', self.on_click)
        self.mpl_connect('button_release_event', self.on_release)
        self.mpl_connect('motion_notify_event', self.on_motion)
        
    def plot_spectrum(self, energies, isotope='Co60'):
        self.ax.clear()
        self.energies = energies
        bin_width = 0.5 
        max_energy = max(energies) if len(energies) > 0 else 2000
        bins = np.arange(0, max_energy + bin_width, bin_width)
        counts, self.bins_edges = np.histogram(energies, bins=bins)
        self.spectrum_data = counts
        self.bins = (self.bins_edges[:-1] + self.bins_edges[1:]) / 2
        
        self.ax.step(self.bins, counts, where='mid', color='b', linewidth=0.5, label='Spectrum')
        self.detect_peaks(isotope)
        
        self.ax.set_xlabel('Energy (keV)')
        self.ax.set_ylabel('Counts')
        self.ax.set_title(f'HPGe Energy Spectrum - {isotope}')
        self.ax.set_yscale('log')
        self.ax.grid(True, alpha=0.3)
        self.draw()
    
    def detect_peaks(self, isotope):
        if isotope not in ISOTOPE_DATABASE: return
        self.detected_peaks = []
        for peak_info in ISOTOPE_DATABASE[isotope]:
            energy = peak_info['energy']
            window_low, window_high = energy - 10.0, energy + 10.0
            mask = (self.bins >= window_low) & (self.bins <= window_high)
            if not np.any(mask): continue
            
            window_counts = self.spectrum_data[mask]
            if len(window_counts) == 0 or np.max(window_counts) < 5: continue
            
            max_idx = np.argmax(window_counts)
            peak_energy = self.bins[mask][max_idx]
            if abs(peak_energy - energy) < 5.0:
                net, gross, bg = calculate_peak_area(self.bins, self.spectrum_data, peak_energy)
                self.detected_peaks.append({
                    'energy': peak_energy, 'counts': window_counts[max_idx],
                    'net_area': net, 'theoretical_energy': energy,
                    'name': peak_info['name'], 'intensity': peak_info['intensity']
                })
                self.ax.axvline(peak_energy, color='r', linestyle='--', alpha=0.5)
                self.ax.text(peak_energy, window_counts[max_idx]*1.5, f'{energy:.1f}', rotation=90, fontsize=8, color='r')

    def on_click(self, event):
        if event.inaxes != self.ax or event.button != 1: return
        self.selecting = True
        self.select_start = event.xdata
        
    def on_motion(self, event):
        if not self.selecting or event.inaxes != self.ax: return
        if self.select_rect: self.select_rect.remove()
        width = abs(event.xdata - self.select_start)
        self.select_rect = Rectangle((min(self.select_start, event.xdata), self.ax.get_ylim()[0]),
                                     width, self.ax.get_ylim()[1], alpha=0.3, color='yellow')
        self.ax.add_patch(self.select_rect)
        self.draw()
    
    def on_release(self, event):
        if not self.selecting: return
        self.selecting = False
        if self.select_rect: 
            self.select_rect.remove()
            self.select_rect = None
            self.draw()
        if event.inaxes != self.ax: return
        
        x_min, x_max = min(self.select_start, event.xdata), max(self.select_start, event.xdata)
        mask = (self.bins >= x_min) & (self.bins <= x_max)
        if np.any(mask):
            total = np.sum(self.spectrum_data[mask])
            centroid = np.sum(self.bins[mask] * self.spectrum_data[mask]) / total if total > 0 else 0
            self.peak_selected.emit(centroid, total)


# ==================================================================
# MARINELLI CONFIG DIALOG
# ==================================================================

class MarinelliConfigDialog(QDialog):
    """Dialog for configuring Marinelli beaker parameters"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Marinelli Beaker Configuration")
        self.resize(500, 400)
        self.marinelli_config = {'enabled': False, 'type': '1000ml', 'fill_material': 'water'}
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        type_group = QGroupBox("Marinelli Beaker Type")
        type_layout = QVBoxLayout()
        self.type_button_group = QButtonGroup()
        self.ml1000_radio = QRadioButton("1000 mL Marinelli Beaker")
        self.ml1000_radio.setChecked(True)
        type_layout.addWidget(self.ml1000_radio)
        type_layout.addWidget(QLabel("  • Outer D: 130mm, Inner d: 85mm, H: 152mm"))
        self.ml450_radio = QRadioButton("450 mL Marinelli Beaker")
        type_layout.addWidget(self.ml450_radio)
        type_layout.addWidget(QLabel("  • Outer D: 114mm, Inner d: 77mm, H: 101mm"))
        self.type_button_group.addButton(self.ml1000_radio, 0)
        self.type_button_group.addButton(self.ml450_radio, 1)
        type_group.setLayout(type_layout)
        layout.addWidget(type_group)
        
        material_group = QGroupBox("Fill Material")
        material_layout = QVBoxLayout()
        self.material_combo = QComboBox()
        self.material_combo.addItems(["water", "soil", "air"])
        material_layout.addWidget(self.material_combo)
        material_group.setLayout(material_layout)
        layout.addWidget(material_group)
        
        button_layout = QHBoxLayout()
        self.ok_button = QPushButton("OK")
        self.ok_button.clicked.connect(self.accept)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        layout.addLayout(button_layout)
        self.setLayout(layout)
    
    def get_configuration(self):
        beaker_type = "1000ml" if self.ml1000_radio.isChecked() else "450ml"
        return {'enabled': True, 'type': beaker_type, 'fill_material': self.material_combo.currentText()}
    
    @staticmethod
    def get_marinelli_config(parent=None):
        dialog = MarinelliConfigDialog(parent)
        result = dialog.exec_()
        return dialog.get_configuration() if result == QDialog.Accepted else None


# ==================================================================
# DISK SOURCE CONFIG DIALOG
# ==================================================================

class DiskSourceConfigDialog(QDialog):
    """Dialog for configuring Disk source parameters"""
    def __init__(self, parent=None, current_config=None):
        super().__init__(parent)
        self.setWindowTitle("📀 Disk Source Configuration")
        self.resize(500, 400)
        self.config = current_config or {
            'enabled': True,
            'radius': 25.0,
            'thickness': 5.0,
            'position_z': -50.0,
            'material': 'water'
        }
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Header
        header = QLabel("<h3>📀 Disk Source Configuration</h3>")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        info_label = QLabel(
            "Configure a flat disk source placed below the detector (Z < 0).\n"
            "The disk is positioned at the specified Z coordinate."
        )
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; margin-bottom: 10px;")
        layout.addWidget(info_label)
        
        # Dimensions group
        dim_group = QGroupBox("Dimensions")
        dim_layout = QGridLayout()
        
        dim_layout.addWidget(QLabel("Radius:"), 0, 0)
        self.radius_spin = QDoubleSpinBox()
        self.radius_spin.setRange(1, 200)
        self.radius_spin.setValue(self.config['radius'])
        self.radius_spin.setSuffix(" mm")
        self.radius_spin.setDecimals(1)
        dim_layout.addWidget(self.radius_spin, 0, 1)
        
        dim_layout.addWidget(QLabel("Thickness:"), 1, 0)
        self.thickness_spin = QDoubleSpinBox()
        self.thickness_spin.setRange(0.1, 50)
        self.thickness_spin.setValue(self.config['thickness'])
        self.thickness_spin.setSuffix(" mm")
        self.thickness_spin.setDecimals(1)
        dim_layout.addWidget(self.thickness_spin, 1, 1)
        
        dim_group.setLayout(dim_layout)
        layout.addWidget(dim_group)
        
        # Position group
        pos_group = QGroupBox("Position")
        pos_layout = QGridLayout()
        
        pos_layout.addWidget(QLabel("Z Position:"), 0, 0)
        self.position_z_spin = QDoubleSpinBox()
        self.position_z_spin.setRange(-500, 0)
        self.position_z_spin.setValue(self.config['position_z'])
        self.position_z_spin.setSuffix(" mm")
        self.position_z_spin.setDecimals(1)
        pos_layout.addWidget(self.position_z_spin, 0, 1)
        
        pos_layout.addWidget(QLabel("(Negative Z = below detector window at Z=0)"), 1, 0, 1, 2)
        
        pos_group.setLayout(pos_layout)
        layout.addWidget(pos_group)
        
        # Material group
        mat_group = QGroupBox("Source Material")
        mat_layout = QVBoxLayout()
        
        self.material_combo = QComboBox()
        self.material_combo.addItems(["water", "soil", "air", "polypropylene"])
        idx = self.material_combo.findText(self.config['material'])
        if idx >= 0:
            self.material_combo.setCurrentIndex(idx)
        mat_layout.addWidget(self.material_combo)
        
        mat_group.setLayout(mat_layout)
        layout.addWidget(mat_group)
        
        # Volume calculation display
        self.volume_label = QLabel("")
        self.volume_label.setStyleSheet("font-weight: bold; color: #0066cc;")
        layout.addWidget(self.volume_label)
        self.update_volume_display()
        
        self.radius_spin.valueChanged.connect(self.update_volume_display)
        self.thickness_spin.valueChanged.connect(self.update_volume_display)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.ok_button = QPushButton("✓ OK")
        self.ok_button.clicked.connect(self.accept)
        self.cancel_button = QPushButton("✗ Cancel")
        self.cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def update_volume_display(self):
        r = self.radius_spin.value() / 10.0  # mm to cm
        h = self.thickness_spin.value() / 10.0  # mm to cm
        vol = 3.14159 * r * r * h
        self.volume_label.setText(f"📐 Calculated Volume: {vol:.2f} cm³")
    
    def get_configuration(self):
        return {
            'enabled': True,
            'radius': self.radius_spin.value(),
            'thickness': self.thickness_spin.value(),
            'position_z': self.position_z_spin.value(),
            'material': self.material_combo.currentText()
        }
    
    @staticmethod
    def get_disk_config(parent=None, current_config=None):
        dialog = DiskSourceConfigDialog(parent, current_config)
        result = dialog.exec_()
        return dialog.get_configuration() if result == QDialog.Accepted else None


# ==================================================================
# CYLINDER SOURCE CONFIG DIALOG
# ==================================================================

class CylinderSourceConfigDialog(QDialog):
    """Dialog for configuring Cylindrical volume source parameters"""
    def __init__(self, parent=None, current_config=None):
        super().__init__(parent)
        self.setWindowTitle("🧪 Cylinder Source Configuration")
        self.resize(550, 550)
        self.config = current_config or {
            'enabled': True,
            'inner_radius': 30.0,
            'outer_radius': 32.0,
            'height': 50.0,
            'wall_thickness': 2.0,
            'bottom_thickness': 2.0,
            'position_z': -60.0,
            'wall_material': 'polypropylene',
            'fill_material': 'water'
        }
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Header
        header = QLabel("<h3>🧪 Cylindrical Volume Source Configuration</h3>")
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        info_label = QLabel(
            "Configure a cylindrical container with walls and sample material.\n"
            "Positioned below the detector (Z < 0). Z position refers to the bottom of the container."
        )
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; margin-bottom: 10px;")
        layout.addWidget(info_label)
        
        # Container dimensions group
        container_group = QGroupBox("Container Dimensions")
        container_layout = QGridLayout()
        
        container_layout.addWidget(QLabel("Outer Radius:"), 0, 0)
        self.outer_radius_spin = QDoubleSpinBox()
        self.outer_radius_spin.setRange(5, 200)
        self.outer_radius_spin.setValue(self.config['outer_radius'])
        self.outer_radius_spin.setSuffix(" mm")
        self.outer_radius_spin.setDecimals(1)
        self.outer_radius_spin.valueChanged.connect(self.update_inner_from_wall)
        container_layout.addWidget(self.outer_radius_spin, 0, 1)
        
        container_layout.addWidget(QLabel("Wall Thickness:"), 1, 0)
        self.wall_thickness_spin = QDoubleSpinBox()
        self.wall_thickness_spin.setRange(0.5, 20)
        self.wall_thickness_spin.setValue(self.config['wall_thickness'])
        self.wall_thickness_spin.setSuffix(" mm")
        self.wall_thickness_spin.setDecimals(1)
        self.wall_thickness_spin.valueChanged.connect(self.update_inner_from_wall)
        container_layout.addWidget(self.wall_thickness_spin, 1, 1)
        
        container_layout.addWidget(QLabel("Inner Radius (calculated):"), 2, 0)
        self.inner_radius_label = QLabel("")
        self.inner_radius_label.setStyleSheet("font-weight: bold;")
        container_layout.addWidget(self.inner_radius_label, 2, 1)
        
        container_layout.addWidget(QLabel("Total Height:"), 3, 0)
        self.height_spin = QDoubleSpinBox()
        self.height_spin.setRange(5, 300)
        self.height_spin.setValue(self.config['height'])
        self.height_spin.setSuffix(" mm")
        self.height_spin.setDecimals(1)
        container_layout.addWidget(self.height_spin, 3, 1)
        
        container_layout.addWidget(QLabel("Bottom Thickness:"), 4, 0)
        self.bottom_thickness_spin = QDoubleSpinBox()
        self.bottom_thickness_spin.setRange(0.5, 20)
        self.bottom_thickness_spin.setValue(self.config['bottom_thickness'])
        self.bottom_thickness_spin.setSuffix(" mm")
        self.bottom_thickness_spin.setDecimals(1)
        container_layout.addWidget(self.bottom_thickness_spin, 4, 1)
        
        container_group.setLayout(container_layout)
        layout.addWidget(container_group)
        
        # Position group
        pos_group = QGroupBox("Position")
        pos_layout = QGridLayout()
        
        pos_layout.addWidget(QLabel("Z Position (bottom):"), 0, 0)
        self.position_z_spin = QDoubleSpinBox()
        self.position_z_spin.setRange(-500, 0)
        self.position_z_spin.setValue(self.config['position_z'])
        self.position_z_spin.setSuffix(" mm")
        self.position_z_spin.setDecimals(1)
        pos_layout.addWidget(self.position_z_spin, 0, 1)
        
        pos_layout.addWidget(QLabel("(Negative Z = below detector window at Z=0)"), 1, 0, 1, 2)
        
        pos_group.setLayout(pos_layout)
        layout.addWidget(pos_group)
        
        # Materials group
        mat_group = QGroupBox("Materials")
        mat_layout = QGridLayout()
        
        mat_layout.addWidget(QLabel("Container Wall Material:"), 0, 0)
        self.wall_material_combo = QComboBox()
        self.wall_material_combo.addItems(["polypropylene", "hdpe", "glass", "pvc", "perspex", "aluminum"])
        idx = self.wall_material_combo.findText(self.config['wall_material'])
        if idx >= 0:
            self.wall_material_combo.setCurrentIndex(idx)
        mat_layout.addWidget(self.wall_material_combo, 0, 1)
        
        mat_layout.addWidget(QLabel("Sample Fill Material:"), 1, 0)
        self.fill_material_combo = QComboBox()
        self.fill_material_combo.addItems(["water", "soil", "air"])
        idx = self.fill_material_combo.findText(self.config['fill_material'])
        if idx >= 0:
            self.fill_material_combo.setCurrentIndex(idx)
        mat_layout.addWidget(self.fill_material_combo, 1, 1)
        
        mat_group.setLayout(mat_layout)
        layout.addWidget(mat_group)
        
        # Volume calculation display
        self.volume_label = QLabel("")
        self.volume_label.setStyleSheet("font-weight: bold; color: #0066cc;")
        layout.addWidget(self.volume_label)
        
        self.update_inner_from_wall()
        
        self.height_spin.valueChanged.connect(self.update_volume_display)
        self.bottom_thickness_spin.valueChanged.connect(self.update_volume_display)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.ok_button = QPushButton("✓ OK")
        self.ok_button.clicked.connect(self.accept)
        self.cancel_button = QPushButton("✗ Cancel")
        self.cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def update_inner_from_wall(self):
        outer_r = self.outer_radius_spin.value()
        wall = self.wall_thickness_spin.value()
        inner_r = outer_r - wall
        self.inner_radius_label.setText(f"{inner_r:.1f} mm")
        self.update_volume_display()
    
    def update_volume_display(self):
        outer_r = self.outer_radius_spin.value() / 10.0  # mm to cm
        wall = self.wall_thickness_spin.value() / 10.0
        inner_r = outer_r - wall
        h = self.height_spin.value() / 10.0
        bottom = self.bottom_thickness_spin.value() / 10.0
        
        sample_h = h - bottom
        sample_vol = 3.14159 * inner_r * inner_r * sample_h
        self.volume_label.setText(f"📐 Sample Volume: {sample_vol:.2f} cm³ (Height: {sample_h*10:.1f} mm)")
    
    def get_configuration(self):
        outer_r = self.outer_radius_spin.value()
        wall = self.wall_thickness_spin.value()
        return {
            'enabled': True,
            'inner_radius': outer_r - wall,
            'outer_radius': outer_r,
            'height': self.height_spin.value(),
            'wall_thickness': wall,
            'bottom_thickness': self.bottom_thickness_spin.value(),
            'position_z': self.position_z_spin.value(),
            'wall_material': self.wall_material_combo.currentText(),
            'fill_material': self.fill_material_combo.currentText()
        }
    
    @staticmethod
    def get_cylinder_config(parent=None, current_config=None):
        dialog = CylinderSourceConfigDialog(parent, current_config)
        result = dialog.exec_()
        return dialog.get_configuration() if result == QDialog.Accepted else None


# ==================================================================
# SIMULATION WORKER
# ==================================================================

class SimulationWorker(QThread):
    progress = pyqtSignal(int)
    finished = pyqtSignal(str)
    output = pyqtSignal(str)
    
    def __init__(self, geant4_path, sim_path, config_file, num_events):
        super().__init__()
        self.geant4_path = geant4_path
        self.sim_path = sim_path
        self.config_file = config_file
        self.num_events = num_events
        self._is_running = True
        
    def run(self):
        macro_path = os.path.join(self.sim_path, 'run_temp.mac')
        with open(macro_path, 'w') as f:
            f.write(f'/run/initialize\n/run/beamOn {self.num_events}\n')
        
        script_path = os.path.join(self.sim_path, 'run_sim.sh')
        with open(script_path, 'w') as f:
            f.write(f'#!/bin/bash\nsource {self.geant4_path}\ncd {self.sim_path}\nexec ./hpge_sim run_temp.mac\n')
        os.chmod(script_path, 0o755)
        
        process = None
        try:
            self.output.emit(f"Running simulation with {self.num_events} events...")
            process = subprocess.Popen(['bash', script_path], stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, universal_newlines=True,
                                     cwd=self.sim_path, preexec_fn=os.setsid)
            self.process = process
            for line in process.stdout:
                if not self._is_running: break
                self.output.emit(line.strip())
            
            if self._is_running:
                process.wait()
                if process.returncode == 0: self.finished.emit("✓ Simulation completed successfully!")
                else: self.finished.emit(f"✗ Simulation failed (Code {process.returncode})")
            else: self.finished.emit("⚠ Stopped by user.")
        except Exception as e:
            self.finished.emit(f"✗ Error: {str(e)}")
        finally:
            if process and process.poll() is None:
                try: os.killpg(os.getpgid(process.pid), subprocess.signal.SIGTERM)
                except: pass

    def stop(self):
        self._is_running = False
        if hasattr(self, 'process') and self.process:
            import signal
            try: os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            except: pass


# ==================================================================
# DETECTOR GEOMETRY WIDGET (Canberra Mirion S/N:21206)
# ==================================================================

class DetectorGeometryWidget(QWidget):
    """Widget for controlling HPGe detector geometry parameters.
    Based on Canberra Mirion S/N:21206 characterization sheet."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
    
    def init_ui(self):
        main_layout = QVBoxLayout()
        
        title = QLabel("<h3>HPGe Detector Geometry - Canberra Mirion S/N: 21206</h3>")
        title.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title)
        
        info = QLabel(
            "Real geometry from characterization sheet. "
            "All dimensions in mm unless noted. Modify to match your detector."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color: #666; margin-bottom: 8px;")
        main_layout.addWidget(info)
        
        # === CRYSTAL GROUP ===
        crystal_group = QGroupBox("Crystal (Ge)")
        crystal_layout = QGridLayout()
        
        crystal_layout.addWidget(QLabel("Crystal Diameter:"), 0, 0)
        self.crystal_diameter = QDoubleSpinBox()
        self.crystal_diameter.setRange(10, 120)
        self.crystal_diameter.setValue(60.5)
        self.crystal_diameter.setSuffix(" mm")
        self.crystal_diameter.setDecimals(1)
        crystal_layout.addWidget(self.crystal_diameter, 0, 1)
        
        crystal_layout.addWidget(QLabel("Crystal Length:"), 0, 2)
        self.crystal_length = QDoubleSpinBox()
        self.crystal_length.setRange(10, 120)
        self.crystal_length.setValue(45.4)
        self.crystal_length.setSuffix(" mm")
        self.crystal_length.setDecimals(1)
        crystal_layout.addWidget(self.crystal_length, 0, 3)
        
        crystal_layout.addWidget(QLabel("Hole Diameter:"), 1, 0)
        self.hole_diameter = QDoubleSpinBox()
        self.hole_diameter.setRange(1, 30)
        self.hole_diameter.setValue(9.5)
        self.hole_diameter.setSuffix(" mm")
        self.hole_diameter.setDecimals(1)
        crystal_layout.addWidget(self.hole_diameter, 1, 1)
        
        crystal_layout.addWidget(QLabel("Hole Depth:"), 1, 2)
        self.hole_depth = QDoubleSpinBox()
        self.hole_depth.setRange(1, 80)
        self.hole_depth.setValue(19.0)
        self.hole_depth.setSuffix(" mm")
        self.hole_depth.setDecimals(1)
        crystal_layout.addWidget(self.hole_depth, 1, 3)
        
        crystal_group.setLayout(crystal_layout)
        main_layout.addWidget(crystal_group)
        
        # === DEAD LAYERS GROUP ===
        dead_group = QGroupBox("Dead Layers / Electrodes")
        dead_layout = QGridLayout()
        
        dead_layout.addWidget(QLabel("Front Dead Layer:"), 0, 0)
        self.front_dead_layer = QDoubleSpinBox()
        self.front_dead_layer.setRange(0, 100)
        self.front_dead_layer.setValue(0.4)
        self.front_dead_layer.setSuffix(" um")
        self.front_dead_layer.setDecimals(2)
        dead_layout.addWidget(self.front_dead_layer, 0, 1)
        
        dead_layout.addWidget(QLabel("Outer Electrode:"), 0, 2)
        self.outer_electrode = QDoubleSpinBox()
        self.outer_electrode.setRange(0, 5)
        self.outer_electrode.setValue(0.5)
        self.outer_electrode.setSuffix(" mm")
        self.outer_electrode.setDecimals(2)
        dead_layout.addWidget(self.outer_electrode, 0, 3)
        
        dead_layout.addWidget(QLabel("Inner Electrode:"), 1, 0)
        self.inner_electrode = QDoubleSpinBox()
        self.inner_electrode.setRange(0, 100)
        self.inner_electrode.setValue(0.3)
        self.inner_electrode.setSuffix(" um eq.Ge")
        self.inner_electrode.setDecimals(2)
        dead_layout.addWidget(self.inner_electrode, 1, 1)
        
        dead_group.setLayout(dead_layout)
        main_layout.addWidget(dead_group)
        
        # === ENDCAP / WINDOW GROUP ===
        endcap_group = QGroupBox("Endcap & Window")
        endcap_layout = QGridLayout()
        
        endcap_layout.addWidget(QLabel("Endcap Diameter:"), 0, 0)
        self.endcap_diameter = QDoubleSpinBox()
        self.endcap_diameter.setRange(30, 150)
        self.endcap_diameter.setValue(76.2)
        self.endcap_diameter.setSuffix(" mm")
        self.endcap_diameter.setDecimals(1)
        endcap_layout.addWidget(self.endcap_diameter, 0, 1)
        
        endcap_layout.addWidget(QLabel("Endcap Material:"), 0, 2)
        self.endcap_material = QComboBox()
        self.endcap_material.addItems(["Al", "Cu", "Steel"])
        endcap_layout.addWidget(self.endcap_material, 0, 3)
        
        endcap_layout.addWidget(QLabel("Endcap Wall:"), 1, 0)
        self.endcap_wall = QDoubleSpinBox()
        self.endcap_wall.setRange(0.1, 10)
        self.endcap_wall.setValue(1.5)
        self.endcap_wall.setSuffix(" mm")
        self.endcap_wall.setDecimals(1)
        endcap_layout.addWidget(self.endcap_wall, 1, 1)
        
        endcap_layout.addWidget(QLabel("Window Thickness:"), 1, 2)
        self.window_thickness = QDoubleSpinBox()
        self.window_thickness.setRange(0.01, 5)
        self.window_thickness.setValue(0.6)
        self.window_thickness.setSuffix(" mm")
        self.window_thickness.setDecimals(2)
        endcap_layout.addWidget(self.window_thickness, 1, 3)
        
        endcap_layout.addWidget(QLabel("Window Material:"), 2, 0)
        self.window_material = QComboBox()
        self.window_material.addItems(["carbon_epoxy", "Al", "Be", "Mylar"])
        endcap_layout.addWidget(self.window_material, 2, 1)
        
        endcap_layout.addWidget(QLabel("Ge to Endcap:"), 2, 2)
        self.ge_to_endcap = QDoubleSpinBox()
        self.ge_to_endcap.setRange(0.5, 30)
        self.ge_to_endcap.setValue(6.0)
        self.ge_to_endcap.setSuffix(" mm")
        self.ge_to_endcap.setDecimals(1)
        endcap_layout.addWidget(self.ge_to_endcap, 2, 3)
        
        endcap_group.setLayout(endcap_layout)
        main_layout.addWidget(endcap_group)
        
        # === HOLDER GROUP ===
        holder_group = QGroupBox("Holder")
        holder_layout = QGridLayout()
        
        holder_layout.addWidget(QLabel("Holder Outer Dia:"), 0, 0)
        self.holder_diameter = QDoubleSpinBox()
        self.holder_diameter.setRange(30, 120)
        self.holder_diameter.setValue(69.0)
        self.holder_diameter.setSuffix(" mm")
        self.holder_diameter.setDecimals(1)
        holder_layout.addWidget(self.holder_diameter, 0, 1)
        
        holder_layout.addWidget(QLabel("Holder Length:"), 0, 2)
        self.holder_length = QDoubleSpinBox()
        self.holder_length.setRange(30, 200)
        self.holder_length.setValue(95.2)
        self.holder_length.setSuffix(" mm")
        self.holder_length.setDecimals(1)
        holder_layout.addWidget(self.holder_length, 0, 3)
        
        holder_layout.addWidget(QLabel("Holder Material:"), 1, 0)
        self.holder_material = QComboBox()
        self.holder_material.addItems(["Al", "Cu", "Steel"])
        holder_layout.addWidget(self.holder_material, 1, 1)
        
        holder_layout.addWidget(QLabel("HDPE Thickness:"), 1, 2)
        self.hdpe_thickness = QDoubleSpinBox()
        self.hdpe_thickness.setRange(0, 20)
        self.hdpe_thickness.setValue(6.5)
        self.hdpe_thickness.setSuffix(" mm")
        self.hdpe_thickness.setDecimals(1)
        holder_layout.addWidget(self.hdpe_thickness, 1, 3)
        
        holder_group.setLayout(holder_layout)
        main_layout.addWidget(holder_group)
        
        # === ACTIVE VOLUME DISPLAY ===
        self.volume_label = QLabel("")
        self.volume_label.setStyleSheet(
            "font-weight: bold; color: #1565C0; padding: 8px; "
            "background-color: #E3F2FD; border-radius: 5px;"
        )
        main_layout.addWidget(self.volume_label)
        
        # Connect signals for live volume calculation
        self.crystal_diameter.valueChanged.connect(self.update_volume_display)
        self.crystal_length.valueChanged.connect(self.update_volume_display)
        self.hole_diameter.valueChanged.connect(self.update_volume_display)
        self.hole_depth.valueChanged.connect(self.update_volume_display)
        self.outer_electrode.valueChanged.connect(self.update_volume_display)
        self.inner_electrode.valueChanged.connect(self.update_volume_display)
        self.front_dead_layer.valueChanged.connect(self.update_volume_display)
        
        # Preset buttons
        preset_layout = QHBoxLayout()
        
        canberra_btn = QPushButton("Reset to Canberra S/N:21206")
        canberra_btn.clicked.connect(self.load_canberra_preset)
        canberra_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; padding: 6px;")
        preset_layout.addWidget(canberra_btn)
        
        fluka_btn = QPushButton("Load FLUKA Defaults")
        fluka_btn.clicked.connect(self.load_fluka_preset)
        fluka_btn.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold; padding: 6px;")
        preset_layout.addWidget(fluka_btn)
        
        main_layout.addLayout(preset_layout)
        
        main_layout.addStretch()
        self.setLayout(main_layout)
        
        self.update_volume_display()
    
    def update_volume_display(self):
        """Calculate and display approximate active volume"""
        import math
        d = self.crystal_diameter.value()
        h = self.crystal_length.value()
        hole_d = self.hole_diameter.value()
        hole_dep = self.hole_depth.value()
        outer_el = self.outer_electrode.value()  # mm
        inner_el = self.inner_electrode.value() / 1000.0  # um -> mm
        front_dl = self.front_dead_layer.value() / 1000.0  # um -> mm
        
        r = d / 2.0 - outer_el
        r_hole = hole_d / 2.0 + inner_el
        
        # Upper part (solid above hole)
        upper_h = h - hole_dep - front_dl
        upper_vol = math.pi * r * r * upper_h / 1000.0  # mm^3 to cm^3
        
        # Lower part (annular around hole)
        lower_h = hole_dep - front_dl
        lower_vol = math.pi * (r * r - r_hole * r_hole) * lower_h / 1000.0
        
        total_vol = upper_vol + lower_vol
        mass = total_vol * 5.323  # Ge density g/cm3
        
        self.volume_label.setText(
            f"Active Volume: {total_vol:.2f} cm3  |  "
            f"Active Mass: {mass:.1f} g  |  "
            f"Relative Efficiency: ~{total_vol / 1.15:.0f}%"
        )
    
    def load_canberra_preset(self):
        """Load Canberra Mirion S/N:21206 values"""
        self.crystal_diameter.setValue(60.5)
        self.crystal_length.setValue(45.4)
        self.hole_diameter.setValue(9.5)
        self.hole_depth.setValue(19.0)
        self.front_dead_layer.setValue(0.4)
        self.outer_electrode.setValue(0.5)
        self.inner_electrode.setValue(0.3)
        self.endcap_diameter.setValue(76.2)
        self.endcap_wall.setValue(1.5)
        self.window_thickness.setValue(0.6)
        self.window_material.setCurrentText("carbon_epoxy")
        self.endcap_material.setCurrentText("Al")
        self.ge_to_endcap.setValue(6.0)
        self.holder_diameter.setValue(69.0)
        self.holder_length.setValue(95.2)
        self.holder_material.setCurrentText("Al")
        self.hdpe_thickness.setValue(6.5)
    
    def load_fluka_preset(self):
        """Load original FLUKA default values"""
        self.crystal_diameter.setValue(65.0)
        self.crystal_length.setValue(45.3)
        self.hole_diameter.setValue(12.0)
        self.hole_depth.setValue(15.0)
        self.front_dead_layer.setValue(900.0)  # 0.9 mm = 900 um
        self.outer_electrode.setValue(1.0)
        self.inner_electrode.setValue(1000.0)  # 1 mm = 1000 um
        self.endcap_diameter.setValue(76.2)
        self.endcap_wall.setValue(1.5)
        self.window_thickness.setValue(0.6)
        self.window_material.setCurrentText("Al")
        self.endcap_material.setCurrentText("Al")
        self.ge_to_endcap.setValue(7.5)
        self.holder_diameter.setValue(69.0)
        self.holder_length.setValue(95.2)
        self.holder_material.setCurrentText("Al")
        self.hdpe_thickness.setValue(0.0)
    
    def get_geometry_config_lines(self):
        """Generate config.txt lines for detector geometry"""
        lines = []
        lines.append(f"# Detector Geometry (Canberra Mirion)")
        lines.append(f"crystal_diameter = {self.crystal_diameter.value()}")
        lines.append(f"crystal_length = {self.crystal_length.value()}")
        lines.append(f"hole_diameter = {self.hole_diameter.value()}")
        lines.append(f"hole_depth = {self.hole_depth.value()}")
        lines.append(f"ge_dead_layer = {self.outer_electrode.value()}")
        lines.append(f"ge_dead_layer_front = {self.front_dead_layer.value() / 1000.0:.6f}")
        lines.append(f"li_dead_layer = {self.inner_electrode.value() / 1000.0:.6f}")
        lines.append(f"li_dead_layer_front = {self.inner_electrode.value() / 1000.0:.6f}")
        lines.append(f"al_window_thickness = {self.window_thickness.value()}")
        lines.append(f"al_cup_thickness = {self.endcap_wall.value()}")
        lines.append(f"vacuum_gap = {self.ge_to_endcap.value()}")
        lines.append(f"endcap_diameter = {self.endcap_diameter.value()}")
        lines.append(f"endcap_material = {self.endcap_material.currentText()}")
        lines.append(f"window_material = {self.window_material.currentText()}")
        lines.append(f"ge_to_endcap_distance = {self.ge_to_endcap.value()}")
        lines.append(f"holder_outer_diameter = {self.holder_diameter.value()}")
        lines.append(f"holder_total_length = {self.holder_length.value()}")
        lines.append(f"holder_material = {self.holder_material.currentText()}")
        lines.append(f"hdpe_thickness = {self.hdpe_thickness.value()}")
        return lines


# ==================================================================
# CONFIGURATION WIDGET
# ==================================================================

class ConfigurationWidget(QWidget):
    """Widget for complete simulation configuration"""
    
    def __init__(self):
        super().__init__()
        self.marinelli_config = None
        self.disk_config = None
        self.cylinder_config = None
        self.geometry_widget = None  # Will be set externally
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        
        mode_group = QGroupBox("Source Configuration Mode")
        mode_layout = QVBoxLayout()
        self.single_isotope_radio = QCheckBox("Single Isotope")
        self.single_isotope_radio.setChecked(True)
        self.single_isotope_radio.toggled.connect(self.on_mode_changed)
        mode_layout.addWidget(self.single_isotope_radio)
        self.multi_isotope_radio = QCheckBox("Multi-Isotope Mixture")
        self.multi_isotope_radio.toggled.connect(self.on_mode_changed)
        mode_layout.addWidget(self.multi_isotope_radio)
        self.custom_energy_radio = QCheckBox("Custom Energy List")
        self.custom_energy_radio.toggled.connect(self.on_mode_changed)
        mode_layout.addWidget(self.custom_energy_radio)
        mode_group.setLayout(mode_layout)
        layout.addWidget(mode_group)
        
        self.single_isotope_widget = self.create_single_isotope_widget()
        layout.addWidget(self.single_isotope_widget)
        self.multi_isotope_widget = self.create_multi_isotope_widget()
        self.multi_isotope_widget.hide()
        layout.addWidget(self.multi_isotope_widget)
        self.custom_energy_widget = self.create_custom_energy_widget()
        self.custom_energy_widget.hide()
        layout.addWidget(self.custom_energy_widget)
        
        geom_group = QGroupBox("Source Geometry")
        geom_layout = QVBoxLayout()
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("Source Type:"))
        self.source_type = QComboBox()
        self.source_type.addItems(['point', 'disk', 'volume', 'marinelli'])
        self.source_type.currentTextChanged.connect(self.on_source_type_changed)
        type_layout.addWidget(self.source_type)
        
        # Configuration buttons for each source type
        self.marinelli_button = QPushButton("⚙️ Configure Marinelli...")
        self.marinelli_button.clicked.connect(self.configure_marinelli)
        self.marinelli_button.setVisible(False)
        type_layout.addWidget(self.marinelli_button)
        
        self.disk_button = QPushButton("⚙️ Configure Disk...")
        self.disk_button.clicked.connect(self.configure_disk)
        self.disk_button.setVisible(False)
        type_layout.addWidget(self.disk_button)
        
        self.cylinder_button = QPushButton("⚙️ Configure Cylinder...")
        self.cylinder_button.clicked.connect(self.configure_cylinder)
        self.cylinder_button.setVisible(False)
        type_layout.addWidget(self.cylinder_button)
        
        geom_layout.addLayout(type_layout)
        
        # Status labels for each source type
        self.marinelli_status = QLabel("")
        self.marinelli_status.setStyleSheet("color: green; font-weight: bold;")
        self.marinelli_status.setVisible(False)
        geom_layout.addWidget(self.marinelli_status)
        
        self.disk_status = QLabel("")
        self.disk_status.setStyleSheet("color: #0066cc; font-weight: bold;")
        self.disk_status.setVisible(False)
        geom_layout.addWidget(self.disk_status)
        
        self.cylinder_status = QLabel("")
        self.cylinder_status.setStyleSheet("color: #9933cc; font-weight: bold;")
        self.cylinder_status.setVisible(False)
        geom_layout.addWidget(self.cylinder_status)
        
        self.position_group = QGroupBox("Position (cm)")
        pos_layout = QGridLayout()
        self.source_x = QDoubleSpinBox()
        self.source_x.setRange(-100, 100)
        self.source_y = QDoubleSpinBox()
        self.source_y.setRange(-100, 100)
        self.source_z = QDoubleSpinBox()
        self.source_z.setRange(-100, 100)
        self.source_z.setValue(-5.0)
        pos_layout.addWidget(QLabel("X:"), 0, 0); pos_layout.addWidget(self.source_x, 0, 1)
        pos_layout.addWidget(QLabel("Y:"), 0, 2); pos_layout.addWidget(self.source_y, 0, 3)
        pos_layout.addWidget(QLabel("Z:"), 1, 0); pos_layout.addWidget(self.source_z, 1, 1)
        self.position_group.setLayout(pos_layout)
        geom_layout.addWidget(self.position_group)
        
        self.dimensions_group = QGroupBox("Dimensions (cm)")
        dim_layout = QGridLayout()
        self.source_radius = QDoubleSpinBox()
        self.source_radius.setRange(0.1, 50)
        self.source_radius.setValue(1.0)
        dim_layout.addWidget(QLabel("Radius:"), 0, 0); dim_layout.addWidget(self.source_radius, 0, 1)
        self.source_height = QDoubleSpinBox()
        self.source_height.setRange(0.1, 50)
        self.source_height.setValue(1.0)
        self.source_height.setEnabled(False)
        dim_layout.addWidget(QLabel("Height:"), 0, 2); dim_layout.addWidget(self.source_height, 0, 3)
        self.dimensions_group.setLayout(dim_layout)
        self.dimensions_group.setVisible(False)
        geom_layout.addWidget(self.dimensions_group)
        
        geom_group.setLayout(geom_layout)
        layout.addWidget(geom_group)
        
        det_group = QGroupBox("Detector Parameters")
        det_layout = QGridLayout()
        self.energy_resolution = QDoubleSpinBox()
        self.energy_resolution.setRange(0.1, 10)
        self.energy_resolution.setValue(0.5)
        self.energy_resolution.setSuffix(" keV")
        det_layout.addWidget(QLabel("Energy Resolution:"), 0, 0)
        det_layout.addWidget(self.energy_resolution, 0, 1)
        
        self.energy_threshold = QDoubleSpinBox()
        self.energy_threshold.setRange(0, 100)
        self.energy_threshold.setValue(10)
        self.energy_threshold.setSuffix(" keV")
        det_layout.addWidget(QLabel("Energy Threshold:"), 0, 2)
        det_layout.addWidget(self.energy_threshold, 0, 3)
        
        self.coinc_window = QDoubleSpinBox()
        self.coinc_window.setRange(0.001, 100)
        self.coinc_window.setValue(1.0)
        self.coinc_window.setSuffix(" µs")
        det_layout.addWidget(QLabel("Coinc. Window:"), 1, 0)
        det_layout.addWidget(self.coinc_window, 1, 1)
        det_group.setLayout(det_layout)
        layout.addWidget(det_group)
        
        tcs_group = QGroupBox("TCS (True Coincidence Summing)")
        tcs_layout = QVBoxLayout()
        self.enable_tcs = QCheckBox("Enable TCS")
        self.enable_tcs.setChecked(True)
        tcs_layout.addWidget(self.enable_tcs)
        
        tcs_window_layout = QHBoxLayout()
        tcs_window_layout.addWidget(QLabel("TCS Window:"))
        self.tcs_window = QDoubleSpinBox()
        self.tcs_window.setRange(0.001, 10)
        self.tcs_window.setValue(0.1)
        self.tcs_window.setSuffix(" µs")
        tcs_window_layout.addWidget(self.tcs_window)
        tcs_layout.addLayout(tcs_window_layout)
        tcs_group.setLayout(tcs_layout)
        layout.addWidget(tcs_group)
        
        fast_group = QGroupBox("Performance")
        fast_layout = QHBoxLayout()
        self.fast_mode = QCheckBox("Fast Mode (reduced physics)")
        fast_layout.addWidget(self.fast_mode)
        fast_group.setLayout(fast_layout)
        layout.addWidget(fast_group)
        
        self.save_btn = QPushButton("💾 Save Configuration")
        self.save_btn.clicked.connect(self.save_config)
        layout.addWidget(self.save_btn)
        
        self.setLayout(layout)

    def create_single_isotope_widget(self):
        widget = QWidget()
        layout = QVBoxLayout()
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Isotope:"))
        self.single_isotope_combo = QComboBox()
        self.single_isotope_combo.addItems(list(ISOTOPE_DATABASE.keys()))
        h_layout.addWidget(self.single_isotope_combo)
        layout.addLayout(h_layout)
        widget.setLayout(layout)
        return widget

    def create_multi_isotope_widget(self):
        widget = QWidget()
        layout = QHBoxLayout()
        
        avail_layout = QVBoxLayout()
        avail_layout.addWidget(QLabel("Available:"))
        self.available_isotopes = QListWidget()
        self.available_isotopes.addItems(list(ISOTOPE_DATABASE.keys()))
        avail_layout.addWidget(self.available_isotopes)
        layout.addLayout(avail_layout)
        
        btn_layout = QVBoxLayout()
        add_btn = QPushButton("→")
        add_btn.clicked.connect(self.add_iso)
        rem_btn = QPushButton("←")
        rem_btn.clicked.connect(self.rem_iso)
        btn_layout.addStretch(); btn_layout.addWidget(add_btn); btn_layout.addWidget(rem_btn); btn_layout.addStretch()
        layout.addLayout(btn_layout)
        
        sel_layout = QVBoxLayout()
        sel_layout.addWidget(QLabel("Selected:"))
        self.selected_isotopes = QListWidget()
        sel_layout.addWidget(self.selected_isotopes)
        layout.addLayout(sel_layout)
        
        widget.setLayout(layout)
        return widget

    def create_custom_energy_widget(self):
        widget = QWidget()
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Custom Energies and Intensities:"))
        self.energy_table = QTableWidget(0, 2)
        self.energy_table.setHorizontalHeaderLabels(["Energy (keV)", "Intensity (%)"])
        self.energy_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        layout.addWidget(self.energy_table)
        
        btn_layout = QHBoxLayout()
        add_row_btn = QPushButton("+ Add Row")
        add_row_btn.clicked.connect(self.add_energy_row)
        del_row_btn = QPushButton("- Delete Row")
        del_row_btn.clicked.connect(self.delete_energy_row)
        btn_layout.addWidget(add_row_btn)
        btn_layout.addWidget(del_row_btn)
        layout.addLayout(btn_layout)
        
        widget.setLayout(layout)
        return widget

    def add_iso(self): 
        if self.available_isotopes.currentItem(): self.selected_isotopes.addItem(self.available_isotopes.currentItem().text())
    def rem_iso(self): 
        self.selected_isotopes.takeItem(self.selected_isotopes.currentRow())
    def add_energy_row(self): 
        r = self.energy_table.rowCount(); self.energy_table.insertRow(r)
    def delete_energy_row(self): 
        self.energy_table.removeRow(self.energy_table.currentRow())
    
    def on_mode_changed(self):
        s = self.sender()
        if not s.isChecked(): return
        self.single_isotope_radio.setChecked(s == self.single_isotope_radio)
        self.multi_isotope_radio.setChecked(s == self.multi_isotope_radio)
        self.custom_energy_radio.setChecked(s == self.custom_energy_radio)
        self.single_isotope_widget.setVisible(s == self.single_isotope_radio)
        self.multi_isotope_widget.setVisible(s == self.multi_isotope_radio)
        self.custom_energy_widget.setVisible(s == self.custom_energy_radio)
        
    def on_source_type_changed(self, t):
        # Hide all source-specific buttons and status labels
        self.marinelli_button.setVisible(False)
        self.marinelli_status.setVisible(False)
        self.disk_button.setVisible(False)
        self.disk_status.setVisible(False)
        self.cylinder_button.setVisible(False)
        self.cylinder_status.setVisible(False)
        
        # Show/hide position group based on source type
        self.position_group.setVisible(t == 'point')
        self.dimensions_group.setVisible(False)  # Hide simple dimensions for all advanced types
        
        if t == 'marinelli':
            self.marinelli_button.setVisible(True)
            if self.marinelli_config:
                self.marinelli_status.setVisible(True)
            else:
                self.configure_marinelli()
        elif t == 'disk':
            self.disk_button.setVisible(True)
            if self.disk_config:
                self.disk_status.setVisible(True)
            else:
                self.configure_disk()
        elif t == 'volume':
            self.cylinder_button.setVisible(True)
            if self.cylinder_config:
                self.cylinder_status.setVisible(True)
            else:
                self.configure_cylinder()
        
    def configure_marinelli(self):
        config = MarinelliConfigDialog.get_marinelli_config(self)
        if config:
            self.marinelli_config = config
            self.marinelli_status.setText(f"✓ Marinelli: {config['type']}, {config['fill_material']}")
            self.marinelli_status.setVisible(True)
    
    def configure_disk(self):
        config = DiskSourceConfigDialog.get_disk_config(self, self.disk_config)
        if config:
            self.disk_config = config
            self.disk_status.setText(
                f"✓ Disk: R={config['radius']:.1f}mm, H={config['thickness']:.1f}mm, "
                f"Z={config['position_z']:.1f}mm, {config['material']}"
            )
            self.disk_status.setVisible(True)
    
    def configure_cylinder(self):
        config = CylinderSourceConfigDialog.get_cylinder_config(self, self.cylinder_config)
        if config:
            self.cylinder_config = config
            self.cylinder_status.setText(
                f"✓ Cylinder: R={config['outer_radius']:.1f}mm, H={config['height']:.1f}mm, "
                f"Z={config['position_z']:.1f}mm\n"
                f"   Wall: {config['wall_material']}, Fill: {config['fill_material']}"
            )
            self.cylinder_status.setVisible(True)

    def set_geometry_widget(self, widget):
        """Set reference to the DetectorGeometryWidget"""
        self.geometry_widget = widget
    
    def generate_config_text(self):
        lines = ["# HPGe Config Generated by GUI"]
        
        # Include detector geometry parameters if geometry widget is connected
        if self.geometry_widget is not None:
            lines.extend(self.geometry_widget.get_geometry_config_lines())
            lines.append("")
        
        if self.single_isotope_radio.isChecked():
            lines.append(f"isotope = {self.single_isotope_combo.currentText()}")
        elif self.multi_isotope_radio.isChecked():
            isos = [self.selected_isotopes.item(i).text() for i in range(self.selected_isotopes.count())]
            ens, ints = [], []
            for iso in isos:
                for p in ISOTOPE_DATABASE.get(iso, []):
                    ens.append(p['energy']); ints.append(p['intensity'])
            lines.append("isotope = Custom")
            lines.append(f"custom_energies = {','.join(map(str, ens))}")
            lines.append(f"custom_intensities = {','.join(map(str, ints))}")
        elif self.custom_energy_radio.isChecked():
            ens, ints = [], []
            for r in range(self.energy_table.rowCount()):
                try:
                    ens.append(float(self.energy_table.item(r,0).text()))
                    ints.append(float(self.energy_table.item(r,1).text()))
                except: pass
            lines.append("isotope = Custom")
            lines.append(f"custom_energies = {','.join(map(str, ens))}")
            lines.append(f"custom_intensities = {','.join(map(str, ints))}")
            
        st = self.source_type.currentText()
        lines.append(f"source_type = {st}")
        if st == 'marinelli':
            if self.marinelli_config:
                lines.append(f"marinelli_type = {self.marinelli_config['type']}")
                lines.append(f"marinelli_fill_material = {self.marinelli_config['fill_material']}")
            else:
                lines.append("marinelli_type = 1000ml\nmarinelli_fill_material = water")
        elif st == 'disk':
            if self.disk_config:
                lines.append(f"disk_radius = {self.disk_config['radius']}")
                lines.append(f"disk_thickness = {self.disk_config['thickness']}")
                lines.append(f"disk_position_z = {self.disk_config['position_z']}")
                lines.append(f"disk_material = {self.disk_config['material']}")
            else:
                lines.append("disk_radius = 25.0")
                lines.append("disk_thickness = 5.0")
                lines.append("disk_position_z = -50.0")
                lines.append("disk_material = water")
        elif st == 'volume':
            if self.cylinder_config:
                lines.append(f"cylinder_inner_radius = {self.cylinder_config['inner_radius']}")
                lines.append(f"cylinder_outer_radius = {self.cylinder_config['outer_radius']}")
                lines.append(f"cylinder_height = {self.cylinder_config['height']}")
                lines.append(f"cylinder_wall_thickness = {self.cylinder_config['wall_thickness']}")
                lines.append(f"cylinder_bottom_thickness = {self.cylinder_config['bottom_thickness']}")
                lines.append(f"cylinder_position_z = {self.cylinder_config['position_z']}")
                lines.append(f"cylinder_wall_material = {self.cylinder_config['wall_material']}")
                lines.append(f"cylinder_fill_material = {self.cylinder_config['fill_material']}")
            else:
                lines.append("cylinder_inner_radius = 28.0")
                lines.append("cylinder_outer_radius = 30.0")
                lines.append("cylinder_height = 50.0")
                lines.append("cylinder_wall_thickness = 2.0")
                lines.append("cylinder_bottom_thickness = 2.0")
                lines.append("cylinder_position_z = -60.0")
                lines.append("cylinder_wall_material = polypropylene")
                lines.append("cylinder_fill_material = water")
        else:
            lines.append(f"source_x = {self.source_x.value()}")
            lines.append(f"source_y = {self.source_y.value()}")
            lines.append(f"source_z = {self.source_z.value()}")
            if st in ['disk', 'volume']: lines.append(f"source_radius = {self.source_radius.value()}")
            if st == 'volume': lines.append(f"source_height = {self.source_height.value()}")
            
        lines.append(f"energy_resolution = {self.energy_resolution.value()}")
        lines.append(f"energy_threshold = {self.energy_threshold.value()}")
        lines.append(f"coincidence_window = {self.coinc_window.value()}")
        
        lines.append(f"enable_tcs = {'true' if self.enable_tcs.isChecked() else 'false'}")
        if self.enable_tcs.isChecked():
            lines.append(f"tcs_coincidence_window = {self.tcs_window.value()}")
        
        lines.append(f"fast_mode = {'true' if self.fast_mode.isChecked() else 'false'}")
            
        return '\n'.join(lines)

    def save_config(self, filename=None):
        txt = self.generate_config_text()
        if not txt: return False
        
        if filename is None:
            filename, _ = QFileDialog.getSaveFileName(self, "Save Config", "config.txt", "Text (*.txt)")
            
        if filename:
            try:
                with open(filename, 'w') as f: f.write(txt)
                if self.sender() == self.save_btn:
                    QMessageBox.information(self, "Success", f"Saved to {filename}")
                return True
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))
                return False
        return False


# ==================================================================
# SIMULATION CONTROL WIDGET (Enhanced with embedded spectrum)
# ==================================================================

class SimulationControlWidget(QWidget):
    """Widget pour le contrôle de la simulation avec Mode A/B Intégré et spectre en temps réel"""
    
    # Signal to notify when an isotope result is ready
    isotope_result_ready = pyqtSignal(str, dict)  # isotope name, results data
    # Signal to clear live results when starting new batch
    clear_live_results_signal = pyqtSignal()
    
    def __init__(self):
        super().__init__()
        self.worker = None
        self.config_widget_ref = None
        self.geant4_path = "/usr/local/geant4/install/bin/geant4.sh" 
        self.sim_path = os.path.join(os.getcwd(), "build")
        
        # TCS State Management
        self.tcs_running = False
        self.tcs_step = 0
        
        # Mode A/B isotope info storage
        self.mode_ab_isotope = None
        self.mode_ab_energies = []
        self.mode_ab_intensities = []
        
        # Multi-isotope batch mode
        self.multi_isotope_mode = False
        self.multi_isotope_list = []
        self.multi_isotope_index = 0
        self.multi_isotope_results = {}
        
        self.init_ui()

    def set_config_widget(self, widget):
        self.config_widget_ref = widget
        
    def init_ui(self):
        layout = QVBoxLayout()
        
        # --- PATHS ---
        paths_group = QGroupBox("Paths Configuration")
        paths_layout = QGridLayout()
        self.g4_path_edit = QLineEdit(self.geant4_path)
        paths_layout.addWidget(QLabel("Geant4:"), 0, 0)
        paths_layout.addWidget(self.g4_path_edit, 0, 1)
        
        self.sim_path_edit = QLineEdit(self.sim_path)
        paths_layout.addWidget(QLabel("Build Dir:"), 1, 0)
        paths_layout.addWidget(self.sim_path_edit, 1, 1)
        paths_group.setLayout(paths_layout)
        layout.addWidget(paths_group)
        
        # --- GENERAL PARAMETERS ---
        params_group = QGroupBox("General Simulation Parameters")
        params_layout = QHBoxLayout()
        self.num_events = QSpinBox()
        self.num_events.setRange(100, 1000000000)
        self.num_events.setValue(100000)
        self.num_events.setSingleStep(10000)
        params_layout.addWidget(QLabel("Number of Events:"))
        params_layout.addWidget(self.num_events)
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # --- GEANT4 VISUALIZATION ---
        vis_run_group = QGroupBox("Geant4 Visualization")
        vis_layout = QHBoxLayout()
        self.run_vis_btn = QPushButton("🔬 Run Geant4 Visualization")
        self.run_vis_btn.clicked.connect(self.run_geant4_visualization)
        self.run_vis_btn.setStyleSheet("background-color: #FCE4EC; font-weight: bold; padding: 5px; color: #880E4F;")
        vis_layout.addWidget(self.run_vis_btn)
        
        self.stop_btn = QPushButton("⏹️ Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop_simulation)
        vis_layout.addWidget(self.stop_btn)
        vis_run_group.setLayout(vis_layout)
        layout.addWidget(vis_run_group)

        # --- MODE A/B RUN (Multi-Isotope Batch Only) ---
        tcs_run_group = QGroupBox("Advanced Analysis: Mode A/B (COI Calculation)")
        tcs_run_group.setStyleSheet("QGroupBox { border: 2px solid #4CAF50; margin-top: 10px; } QGroupBox::title { color: #2E7D32; font-weight: bold; }")
        tcs_layout = QVBoxLayout()
        
        tcs_info = QLabel("Mode A: Realistic cascade emission | Mode B: Independent single-gamma emission\n"
                         "Results will be displayed in real-time and saved automatically.")
        tcs_info.setWordWrap(True)
        tcs_info.setStyleSheet("font-style: italic; color: #555; padding: 5px;")
        tcs_layout.addWidget(tcs_info)
        
        # Main batch button (now the only run button)
        self.run_multi_ab_btn = QPushButton("▶️ Run Multi-Isotope Batch Analysis")
        self.run_multi_ab_btn.clicked.connect(self.run_multi_isotope_ab)
        self.run_multi_ab_btn.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50; 
                font-weight: bold; 
                padding: 10px; 
                color: white;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:disabled {
                background-color: #cccccc;
                color: #666666;
            }
        """)
        tcs_layout.addWidget(self.run_multi_ab_btn)
        
        # === EMBEDDED LIVE SPECTRUM (Apex-style) ===
        self.embedded_spectrum = EmbeddedLiveSpectrum(self)
        tcs_layout.addWidget(self.embedded_spectrum)
        
        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setStyleSheet("background-color: #ccc;")
        tcs_layout.addWidget(separator)
        
        # Results management buttons
        results_layout = QHBoxLayout()
        
        self.view_results_btn = QPushButton("📊 View Last Batch Results")
        self.view_results_btn.clicked.connect(self.reopen_batch_results)
        self.view_results_btn.setStyleSheet("background-color: #FFF3E0; font-weight: bold; padding: 5px; color: #E65100;")
        results_layout.addWidget(self.view_results_btn)
        
        self.load_results_btn = QPushButton("📂 Load Results from File")
        self.load_results_btn.clicked.connect(self.load_results_from_file)
        self.load_results_btn.setStyleSheet("background-color: #E8EAF6; font-weight: bold; padding: 5px; color: #3F51B5;")
        results_layout.addWidget(self.load_results_btn)
        
        tcs_layout.addLayout(results_layout)
        
        tcs_run_group.setLayout(tcs_layout)
        layout.addWidget(tcs_run_group)
        
        # --- LOG & PROGRESS ---
        self.progress_bar = QProgressBar()
        layout.addWidget(self.progress_bar)
        
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(150)
        layout.addWidget(self.log_text)
        
        self.setLayout(layout)
    
    def run_geant4_visualization(self):
        """Lance Geant4 avec visualisation OpenGL en utilisant l'environnement correct"""
        if not self.check_paths(): return
        
        sim_path = self.sim_path_edit.text()
        g4_env_path = self.g4_path_edit.text() # Chemin vers geant4.sh
        
        # IMPORTANT: Sauvegarder la configuration avant de lancer la visualisation
        if self.config_widget_ref:
            config_path = os.path.join(sim_path, 'config.txt')
            saved = self.config_widget_ref.save_config(config_path)
            if saved:
                self.log_text.append(f"✓ Configuration saved to: {config_path}")
            else:
                self.log_text.append("⚠ Warning: Could not save configuration")
        
        # 1. Protection du fichier vis.mac existant
        # Si vous voulez utiliser le vis.mac qui est DÉJÀ dans le dossier (celui qui marche en terminal),
        # il ne faut pas le réécrire ici. J'ai commenté cette partie.
        vis_mac_path = os.path.join(sim_path, "vis.mac")
        
        # Décommentez les lignes suivantes SEULEMENT si vous voulez que le GUI génère son propre vis.mac
        """
        vis_content = "# Visualization macro for Geant4\n/run/initialize\n/vis/open OGL 800x600-0+0\n/vis/drawVolume\n/vis/viewer/set/viewpointThetaPhi 45 45\n/vis/viewer/set/autoRefresh true\n/vis/scene/add/trajectories smooth\n/vis/scene/endOfEventAction accumulate\n"
        try:
            with open(vis_mac_path, 'w') as f:
                f.write(vis_content)
            self.log_text.append(f"Created temp visualization macro: {vis_mac_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create visualization macro: {e}")
            return
        """
        
        self.log_text.append("=" * 10)
        self.log_text.append("🔬 Starting Geant4 Visualization...")
        
        try:
            # 2. Construction de la commande corrigée
            # On source geant4.sh D'ABORD, puis on va dans le dossier, puis on lance l'exe sans arguments
            # (puisque votre main.cc lit vis.mac par défaut si aucun argument n'est donné)
            
            cmd = f"source {g4_env_path} && cd {sim_path} && ./hpge_sim"
            
            self.log_text.append(f"Executing: {cmd}")
            
            # Utilisation de Popen avec executable='/bin/bash' pour bien gérer le 'source'
            self.vis_process = subprocess.Popen(
                cmd,
                shell=True,
                executable='/bin/bash',  # Important pour que 'source' fonctionne
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=sim_path
            )
            
            self.log_text.append(f"Geant4 visualization started (PID: {self.vis_process.pid})")
            
            self.run_vis_btn.setEnabled(False)
            self.stop_btn.setEnabled(True)
            
            self.vis_timer = QTimer()
            self.vis_timer.timeout.connect(self.check_visualization_status)
            self.vis_timer.start(1000)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start visualization: {e}")
            self.log_text.append(f"Error: {e}")

    def check_visualization_status(self):
        """Check if visualization process is still running"""
        if hasattr(self, 'vis_process') and self.vis_process:
            poll = self.vis_process.poll()
            if poll is not None:
                self.vis_timer.stop()
                self.run_vis_btn.setEnabled(True)
                self.stop_btn.setEnabled(False)
                self.log_text.append("Geant4 visualization session ended.")

    def run_multi_isotope_ab(self):
        """Lance la séquence Mode A/B pour plusieurs isotopes avec affichage en temps réel"""
        if not self.config_widget_ref:
            QMessageBox.critical(self, "Error", "Configuration widget not linked.")
            return
        if not self.check_paths(): return
        
        # Create isotope selection dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Select Isotopes for Batch Analysis")
        dialog.setMinimumWidth(500)
        dlg_layout = QVBoxLayout()
        
        dlg_layout.addWidget(QLabel("<b>Select isotopes to analyze:</b><br>"
                                    "<i>Each isotope will run Mode A + Mode B simulations</i>"))
        
        scroll = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout()
        
        checkboxes = {}
        for iso in ISOTOPE_DATABASE.keys():
            peaks = ISOTOPE_DATABASE[iso]
            energies = ", ".join([f"{p['energy']:.1f}" for p in peaks])
            cascade_marker = " ⚡" if iso in ['Co60', 'Y88', 'Na22', 'Eu152', 'Ba133'] else ""
            cb = QCheckBox(f"{iso}{cascade_marker} ({energies} keV)")
            checkboxes[iso] = cb
            scroll_layout.addWidget(cb)
        
        scroll_widget.setLayout(scroll_layout)
        scroll.setWidget(scroll_widget)
        scroll.setWidgetResizable(True)
        scroll.setMaximumHeight(300)
        dlg_layout.addWidget(scroll)
        
        btn_layout = QHBoxLayout()
        select_all_btn = QPushButton("Select All")
        select_all_btn.clicked.connect(lambda: [cb.setChecked(True) for cb in checkboxes.values()])
        deselect_all_btn = QPushButton("Deselect All")
        deselect_all_btn.clicked.connect(lambda: [cb.setChecked(False) for cb in checkboxes.values()])
        btn_layout.addWidget(select_all_btn)
        btn_layout.addWidget(deselect_all_btn)
        dlg_layout.addLayout(btn_layout)
        
        info_label = QLabel(f"<i>⚡ = Cascade isotope (TCS effects expected)<br>"
                           f"Events per simulation: {self.num_events.value():,}</i>")
        dlg_layout.addWidget(info_label)
        
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        dlg_layout.addWidget(button_box)
        
        dialog.setLayout(dlg_layout)
        
        if dialog.exec_() != QDialog.Accepted:
            return
        
        selected = [iso for iso, cb in checkboxes.items() if cb.isChecked()]
        
        if not selected:
            QMessageBox.warning(self, "Warning", "No isotopes selected.")
            return
        
        total_sims = len(selected) * 2
        reply = QMessageBox.question(self, "Confirm Batch Run",
            f"You selected {len(selected)} isotopes:\n"
            f"{', '.join(selected)}\n\n"
            f"This will run {total_sims} simulations total.\n"
            f"Live spectrum will be displayed during acquisition.\n\n"
            f"Continue?",
            QMessageBox.Yes | QMessageBox.No)
        
        if reply == QMessageBox.No:
            return
        
        # Initialize multi-isotope batch
        self.multi_isotope_mode = True
        self.multi_isotope_list = selected
        self.multi_isotope_index = 0
        self.multi_isotope_results = {}
        
        self.log_text.clear()
        self.log_text.append("=" * 60)
        self.log_text.append("MULTI-ISOTOPE BATCH MODE A/B ANALYSIS")
        self.log_text.append(f"Isotopes: {', '.join(selected)}")
        self.log_text.append(f"Total simulations: {total_sims}")
        self.log_text.append("=" * 60)
        
        # Clear previous live results if results widget is connected
        self.clear_live_results_signal.emit()
        
        # Show the embedded spectrum
        self.embedded_spectrum.setVisible(True)
        
        # Start first isotope
        self.start_next_isotope_ab()

    def start_next_isotope_ab(self):
        """Start Mode A/B for the next isotope in the batch"""
        if self.multi_isotope_index >= len(self.multi_isotope_list):
            self.finish_multi_isotope_batch()
            return
        
        isotope = self.multi_isotope_list[self.multi_isotope_index]
        total = len(self.multi_isotope_list)
        
        self.log_text.append(f"\n{'='*40}")
        self.log_text.append(f"[{self.multi_isotope_index + 1}/{total}] Starting {isotope}...")
        self.log_text.append(f"{'='*40}")
        
        # Set up for this isotope
        self.config_widget_ref.single_isotope_combo.setCurrentText(isotope)
        self.config_widget_ref.single_isotope_radio.setChecked(True)
        
        isotope_peaks = ISOTOPE_DATABASE[isotope]
        self.mode_ab_isotope = isotope
        self.mode_ab_energies = [p['energy'] for p in isotope_peaks]
        self.mode_ab_intensities = [p['intensity'] for p in isotope_peaks]
        
        # Start Mode A
        self.tcs_running = True
        self.tcs_step = 1
        
        # Update embedded spectrum display
        self.embedded_spectrum.set_current_simulation(isotope, "A", self.multi_isotope_index + 1, total)
        self.embedded_spectrum.clear_spectrum()
        
        self.config_widget_ref.enable_tcs.setChecked(True)
        config_path = os.path.join(self.sim_path_edit.text(), 'config.txt')
        self.config_widget_ref.save_config(config_path)
        
        self.log_text.append(f"Running Mode A for {isotope}...")
        self.start_worker_process(f"Mode A ({isotope})...")

    def finish_multi_isotope_batch(self):
        """Finish the batch and show results"""
        self.multi_isotope_mode = False
        self.log_text.append("\n" + "=" * 60)
        self.log_text.append("✅ BATCH COMPLETE!")
        self.log_text.append("=" * 60)
        
        # Update embedded spectrum to show completion
        self.embedded_spectrum.finish_batch()
        
        # Show summary results
        self.show_multi_isotope_results()
        
        self.reset_tcs_state()
        self.progress_bar.setValue(100)

    def show_multi_isotope_results(self):
        """Display results table for all isotopes in a professional dialog"""
        self.batch_results_dialog = BatchResultsDialog(self, self.multi_isotope_results, 
                                                        self.num_events.value())
        self.batch_results_dialog.show()

    def reopen_batch_results(self):
        """Reopen the batch results dialog - from memory or from saved file"""
        if hasattr(self, 'multi_isotope_results') and self.multi_isotope_results:
            self.show_multi_isotope_results()
            return
        
        save_path = os.path.join(os.getcwd(), 'batch_results_cache.json')
        if os.path.exists(save_path):
            try:
                with open(save_path, 'r') as f:
                    saved_data = json.load(f)
                
                self.multi_isotope_results = saved_data['results']
                beamon = saved_data['beamon']
                timestamp = saved_data.get('timestamp', 'Unknown')
                
                reply = QMessageBox.question(self, "Load Saved Results",
                    f"Found saved results from: {timestamp}\n\n"
                    f"Load these results?",
                    QMessageBox.Yes | QMessageBox.No)
                
                if reply == QMessageBox.Yes:
                    self.batch_results_dialog = BatchResultsDialog(self, self.multi_isotope_results, beamon)
                    self.batch_results_dialog.show()
                    return
            except Exception as e:
                print(f"Error loading saved results: {e}")
        
        QMessageBox.warning(self, "No Results", 
                           "No batch results available.\n"
                           "Run a multi-isotope batch first or check if 'batch_results_cache.json' exists.")

    def load_results_from_file(self):
        """Load batch results from a JSON file"""
        filename, _ = QFileDialog.getOpenFileName(self, "Load Batch Results", 
                                                   "", "JSON Files (*.json);;All Files (*)")
        if not filename:
            return
        
        try:
            with open(filename, 'r') as f:
                saved_data = json.load(f)
            
            if 'results' not in saved_data or 'beamon' not in saved_data:
                QMessageBox.critical(self, "Invalid File", 
                                    "The selected file does not contain valid batch results.")
                return
            
            self.multi_isotope_results = saved_data['results']
            beamon = saved_data['beamon']
            timestamp = saved_data.get('timestamp', 'Unknown')
            
            self.batch_results_dialog = BatchResultsDialog(self, self.multi_isotope_results, beamon)
            self.batch_results_dialog.show()
            
            self.log_text.append(f"✓ Loaded batch results from: {filename}")
            self.log_text.append(f"  Timestamp: {timestamp}")
            self.log_text.append(f"  Isotopes: {len(self.multi_isotope_results)}")
            
        except Exception as e:
            QMessageBox.critical(self, "Load Error", f"Error loading file:\n{str(e)}")

    def start_worker_process(self, start_msg):
        self.run_vis_btn.setEnabled(False)
        self.run_multi_ab_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.progress_bar.setValue(0)
        self.log_text.append(start_msg)
        
        # Start embedded live spectrum
        output_file = os.path.join(self.sim_path_edit.text(), "output.dat")
        self.embedded_spectrum.start_live_update(output_file)
        
        self.worker = SimulationWorker(self.g4_path_edit.text(), self.sim_path_edit.text(), 
                                      'config.txt', self.num_events.value())
        self.worker.finished.connect(self.on_simulation_finished)
        self.worker.output.connect(self.update_log)
        self.worker.start()
        
    def on_simulation_finished(self, message):
        is_success = "successfully" in message
        self.log_text.append(message)
        
        # Stop embedded live spectrum
        self.embedded_spectrum.stop_live_update()
        
        if not self.tcs_running:
            self.stop_btn.setEnabled(False)
            self.run_vis_btn.setEnabled(True)
            self.run_multi_ab_btn.setEnabled(True)
            self.progress_bar.setValue(100 if is_success else 0)
            return

        if not is_success:
            self.log_text.append("SEQUENCE ABORTED DUE TO ERROR.")
            self.reset_tcs_state()
            return
            
        sim_path = self.sim_path_edit.text()
        default_output = os.path.join(sim_path, "output.dat") 
        
        if self.tcs_step == 1:
            # Finished Mode A
            self.log_text.append("Step 1 Finished. Saving output to 'output_mode_a_tcs.txt'...")
            target_a = os.path.join(sim_path, "output_mode_a_tcs.txt")
            try:
                if os.path.exists(target_a): os.remove(target_a)
                shutil.move(default_output, target_a)
            except Exception as e:
                self.log_text.append(f"Error moving file: {e}")
                self.reset_tcs_state()
                return

            # Start Step 2 - Mode B
            self.tcs_step = 2
            energies_str = ",".join([f"{e:.2f}" for e in self.mode_ab_energies])
            intensities_str = ",".join([f"{i:.2f}" for i in self.mode_ab_intensities])
            
            self.log_text.append(f"\nStep 2/2: Configuring Mode B (Independent emission)...")
            
            # Update embedded spectrum for Mode B
            self.embedded_spectrum.set_current_simulation(
                self.mode_ab_isotope, "B", 
                self.multi_isotope_index + 1, len(self.multi_isotope_list)
            )
            self.embedded_spectrum.clear_spectrum()
            
            config_path = os.path.join(sim_path, 'config.txt')
            self.save_mode_b_config(config_path)
            
            self.start_worker_process("Running Mode B (Independent single-gamma emission)...")
            
        elif self.tcs_step == 2:
            # Finished Mode B
            self.log_text.append("Step 2 Finished. Saving output to 'output_mode_b_notcs.txt'...")
            target_b = os.path.join(sim_path, "output_mode_b_notcs.txt")
            try:
                if os.path.exists(target_b): os.remove(target_b)
                shutil.move(default_output, target_b)
            except Exception as e:
                self.log_text.append(f"Error moving file: {e}")
            
            if self.multi_isotope_mode:
                self.store_current_isotope_results(sim_path)
                
                # Emit signal that isotope results are ready (for creating new tab)
                self.isotope_result_ready.emit(self.mode_ab_isotope, 
                                               self.multi_isotope_results.get(self.mode_ab_isotope, {}))
                
                # Move to next isotope
                self.multi_isotope_index += 1
                progress = int((self.multi_isotope_index / len(self.multi_isotope_list)) * 100)
                self.progress_bar.setValue(progress)
                
                # Reset for next isotope
                self.tcs_step = 0
                self.tcs_running = False
                
                # Start next isotope
                self.start_next_isotope_ab()
            else:
                self.config_widget_ref.single_isotope_radio.setChecked(True)
                
                self.log_text.append("\n--- MODE A/B SEQUENCE COMPLETED SUCCESSFULLY ---")
                self.reset_tcs_state()
                self.progress_bar.setValue(100)

    def store_current_isotope_results(self, sim_path):
        """Store the results for the current isotope in multi-isotope mode"""
        isotope = self.mode_ab_isotope
        
        file_a = os.path.join(sim_path, "output_mode_a_tcs.txt")
        file_b = os.path.join(sim_path, "output_mode_b_notcs.txt")
        
        try:
            energies_a = self.load_raw_energies(file_a)
            energies_b = self.load_raw_energies(file_b)
        except Exception as e:
            self.log_text.append(f"Error reading results for {isotope}: {e}")
            return
        
        isotope_peaks = ISOTOPE_DATABASE.get(isotope, [])
        energy_results = []
        
        for i, peak in enumerate(isotope_peaks):
            energy = peak['energy']
            intensity = peak['intensity']
            
            cnt_a = self.get_peak_counts_local(energies_a, energy)
            cnt_b = self.get_peak_counts_local(energies_b, energy)
            
            energy_results.append({
                'energy': energy,
                'intensity': intensity,
                'cnt_a': int(cnt_a),
                'cnt_b': int(cnt_b)
            })
        
        self.multi_isotope_results[isotope] = {
            'energies': energy_results
        }
        
        try:
            archive_a = os.path.join(sim_path, f"output_{isotope}_mode_a.txt")
            archive_b = os.path.join(sim_path, f"output_{isotope}_mode_b.txt")
            shutil.copy(file_a, archive_a)
            shutil.copy(file_b, archive_b)
        except Exception as e:
            self.log_text.append(f"Warning: Could not archive files for {isotope}: {e}")
        
        self.log_text.append(f"✓ Results stored for {isotope}")

    def load_raw_energies(self, path):
        """Load raw energies from output file"""
        data = np.loadtxt(path, comments='#')
        energies = []
        if len(data.shape) < 2: return np.array(energies)
        
        for row in data:
            if len(row) < 3: continue
            n_det = int(row[2])
            for i in range(n_det):
                idx = 4 + i*3
                if idx < len(row): energies.append(row[idx])
        return np.array(energies)

    def get_peak_counts_local(self, energies, peak_energy):
        """Get counts in 3% window around peak energy"""
        w_l = peak_energy * 0.97
        w_h = peak_energy * 1.03
        counts = np.sum((energies >= w_l) & (energies <= w_h))
        return counts

    def save_mode_b_config(self, filepath):
        """Generate and save config for Mode B with Custom mode (independent emission)"""
        cfg = self.config_widget_ref
        
        energies_str = ",".join([f"{e:.2f}" for e in self.mode_ab_energies])
        intensities_str = ",".join([f"{i:.2f}" for i in self.mode_ab_intensities])
        
        lines = ["# HPGe Config Generated by GUI - Mode B (Independent emission)"]
        lines.append("isotope = Custom")
        lines.append(f"custom_energies = {energies_str}")
        lines.append(f"custom_intensities = {intensities_str}")
        
        st = cfg.source_type.currentText()
        lines.append(f"source_type = {st}")
        
        if st == 'marinelli':
            if cfg.marinelli_config:
                lines.append(f"marinelli_type = {cfg.marinelli_config['type']}")
                lines.append(f"marinelli_fill_material = {cfg.marinelli_config['fill_material']}")
            else:
                lines.append("marinelli_type = 1000ml")
                lines.append("marinelli_fill_material = water")
        elif st == 'disk':
            # Use disk_config if available
            if cfg.disk_config:
                lines.append(f"disk_radius = {cfg.disk_config['radius']}")
                lines.append(f"disk_thickness = {cfg.disk_config['thickness']}")
                lines.append(f"disk_position_z = {cfg.disk_config['position_z']}")
                lines.append(f"disk_material = {cfg.disk_config['material']}")
            else:
                # Default values
                lines.append("disk_radius = 25.0")
                lines.append("disk_thickness = 5.0")
                lines.append("disk_position_z = -50.0")
                lines.append("disk_material = water")
        elif st == 'volume':
            # Use cylinder_config if available
            if cfg.cylinder_config:
                lines.append(f"cylinder_inner_radius = {cfg.cylinder_config['inner_radius']}")
                lines.append(f"cylinder_outer_radius = {cfg.cylinder_config['outer_radius']}")
                lines.append(f"cylinder_height = {cfg.cylinder_config['height']}")
                lines.append(f"cylinder_wall_thickness = {cfg.cylinder_config['wall_thickness']}")
                lines.append(f"cylinder_bottom_thickness = {cfg.cylinder_config['bottom_thickness']}")
                lines.append(f"cylinder_position_z = {cfg.cylinder_config['position_z']}")
                lines.append(f"cylinder_wall_material = {cfg.cylinder_config['wall_material']}")
                lines.append(f"cylinder_fill_material = {cfg.cylinder_config['fill_material']}")
            else:
                # Default values
                lines.append("cylinder_inner_radius = 28.0")
                lines.append("cylinder_outer_radius = 30.0")
                lines.append("cylinder_height = 50.0")
                lines.append("cylinder_wall_thickness = 2.0")
                lines.append("cylinder_bottom_thickness = 2.0")
                lines.append("cylinder_position_z = -60.0")
                lines.append("cylinder_wall_material = polypropylene")
                lines.append("cylinder_fill_material = water")
        else:
            # Point source - use position values
            lines.append(f"source_x = {cfg.source_x.value()}")
            lines.append(f"source_y = {cfg.source_y.value()}")
            lines.append(f"source_z = {cfg.source_z.value()}")
        
        lines.append(f"energy_resolution = {cfg.energy_resolution.value()}")
        lines.append(f"energy_threshold = {cfg.energy_threshold.value()}")
        lines.append(f"coincidence_window = {cfg.coinc_window.value()}")
        
        lines.append("enable_tcs = false")
        lines.append(f"fast_mode = {'true' if cfg.fast_mode.isChecked() else 'false'}")
        
        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))

    def reset_tcs_state(self):
        self.tcs_running = False
        self.tcs_step = 0
        self.multi_isotope_mode = False
        self.run_vis_btn.setEnabled(True)
        self.run_multi_ab_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def update_log(self, msg):
        self.log_text.append(msg)
        if "Run terminated." in msg: self.progress_bar.setValue(90)
        
    def stop_simulation(self):
        if self.worker: self.worker.stop()
        if hasattr(self, 'vis_process') and self.vis_process:
            try:
                self.vis_process.terminate()
                self.vis_process = None
                self.log_text.append("Visualization stopped.")
            except:
                pass
        if hasattr(self, 'vis_timer'):
            self.vis_timer.stop()
        
        self.embedded_spectrum.stop_live_update()
        self.run_vis_btn.setEnabled(True)
        self.reset_tcs_state()

    def check_paths(self):
        if not os.path.exists(self.sim_path_edit.text()):
            QMessageBox.critical(self, "Error", "Build directory not found!")
            return False
        return True


# ==================================================================
# ENHANCED RESULTS WIDGET
# ==================================================================

class EnhancedResultsWidget(QWidget):
    """Widget amélioré pour l'analyse interactive des résultats"""
    
    def __init__(self):
        super().__init__()
        self.config_widget_ref = None
        self.sim_widget_ref = None
        self.tab_widget_ref = None
        self.current_isotope = 'Co60'
        self.energies = None
        
        # Storage for live batch results
        self.live_batch_results = {}  # {isotope: {energies: [...]}}
        self.live_results_tab = None
        self.live_results_table = None
        
        self.init_ui()

    def set_references(self, config_widget, sim_widget, tab_widget):
        self.config_widget_ref = config_widget
        self.sim_widget_ref = sim_widget
        self.tab_widget_ref = tab_widget
        
        # Connect to isotope result signal
        if sim_widget:
            sim_widget.isotope_result_ready.connect(self.on_isotope_result_ready)
            sim_widget.clear_live_results_signal.connect(self.clear_live_results)
        
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Toolbar
        toolbar = QHBoxLayout()
        load_btn = QPushButton("📂 Load Last Result")
        load_btn.clicked.connect(self.load_results_dialog)
        toolbar.addWidget(load_btn)
        
        toolbar.addWidget(QLabel("BeamOn:"))
        self.beamon_input = QSpinBox()
        self.beamon_input.setRange(1, int(1e9))
        self.beamon_input.setValue(100000)
        toolbar.addWidget(self.beamon_input)
        
        self.isotope_combo = QComboBox()
        self.isotope_combo.addItems(list(ISOTOPE_DATABASE.keys()) + ['Custom'])
        self.isotope_combo.currentTextChanged.connect(self.on_isotope_changed)
        toolbar.addWidget(QLabel("Iso:")); toolbar.addWidget(self.isotope_combo)
        layout.addLayout(toolbar)
        
        # TCS Analysis Group
        tcs_group = QGroupBox("TCS (True Coincidence Summing) Analysis Results")
        tcs_group.setStyleSheet("QGroupBox { border: 1px solid #28A745; margin-top: 10px; } QGroupBox::title { color: #28A745; }")
        tcs_layout = QVBoxLayout()
        
        btn_row = QHBoxLayout()
        
        self.coi_btn = QPushButton("📊 Load Single Isotope A/B Files")
        self.coi_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.coi_btn.setToolTip("Load output_mode_a_tcs.txt and output_mode_b_notcs.txt")
        self.coi_btn.clicked.connect(self.calculate_and_show_coi)
        btn_row.addWidget(self.coi_btn)
        
        self.multi_coi_btn = QPushButton("📂 Load Multi-Isotope A/B Files")
        self.multi_coi_btn.setStyleSheet("font-weight: bold; padding: 6px; background-color: #E3F2FD;")
        self.multi_coi_btn.setToolTip("Load all output_ISOTOPE_mode_a/b.txt files from batch run")
        self.multi_coi_btn.clicked.connect(self.load_multi_isotope_files)
        btn_row.addWidget(self.multi_coi_btn)
        
        tcs_layout.addLayout(btn_row)
        
        tcs_info = QLabel("Single: 'output_mode_a_tcs.txt' / 'output_mode_b_notcs.txt' | Multi: 'output_ISOTOPE_mode_a/b.txt'")
        tcs_info.setStyleSheet("color: #666; font-size: 10px;")
        tcs_layout.addWidget(tcs_info)
        tcs_group.setLayout(tcs_layout)
        layout.addWidget(tcs_group)

        # Plot Area
        self.canvas = InteractiveSpectrum()
        self.canvas.peak_selected.connect(self.on_peak_selected)
        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(self.canvas)
        
        # Table
        self.peaks_table = QTableWidget(0, 7)
        self.peaks_table.setHorizontalHeaderLabels(["Name", "E_theo", "E_meas", "Counts", "Area(3%)", "Int(%)", "Eff"])
        self.peaks_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        layout.addWidget(self.peaks_table)
        
        # Stats
        self.stats_text = QTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setMaximumHeight(150)
        layout.addWidget(self.stats_text)
        
        self.setLayout(layout)

    def on_isotope_result_ready(self, isotope, results):
        """Called when an isotope simulation is complete - updates the single results tab"""
        # Store results
        self.live_batch_results[isotope] = results
        
        # Update or create the single results tab
        self.update_live_results_tab()

    def update_live_results_tab(self):
        """Update the single live results tab with all completed isotopes"""
        if not self.tab_widget_ref:
            return
        
        tab_name = "📊 Live Batch Results"
        
        # Find existing tab or create new one
        tab_index = -1
        for i in range(self.tab_widget_ref.count()):
            if self.tab_widget_ref.tabText(i) == tab_name:
                tab_index = i
                break
        
        # Get beamon value
        beamon = self.sim_widget_ref.num_events.value() if self.sim_widget_ref else 100000
        cascade_isotopes = ['Co60', 'Y88', 'Na22', 'Eu152', 'Ba133']
        
        if tab_index == -1:
            # Create new tab
            self.live_results_tab = QWidget()
            tab_layout = QVBoxLayout()
            
            # Header
            self.live_header = QLabel("<h2>📊 Live Batch Results</h2>")
            self.live_header.setAlignment(Qt.AlignCenter)
            tab_layout.addWidget(self.live_header)
            
            # Info bar
            self.live_info = QLabel()
            self.live_info.setStyleSheet("background-color: #E3F2FD; padding: 10px; border-radius: 5px;")
            self.live_info.setAlignment(Qt.AlignCenter)
            tab_layout.addWidget(self.live_info)
            
            # Create table
            self.live_results_table = QTableWidget()
            self.live_results_table.setColumnCount(9)
            self.live_results_table.setHorizontalHeaderLabels([
                "Isotope", "Energy\n(keV)", "Intensity\n(%)", 
                "Counts A", "Counts B", "Eff A", "Eff B", "COI", "COI ± σ"
            ])
            self.live_results_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.live_results_table.setAlternatingRowColors(True)
            self.live_results_table.setStyleSheet(
                "QTableWidget { gridline-color: #d0d0d0; }"
                "QTableWidget::item:alternate { background-color: #f5f5f5; }"
            )
            tab_layout.addWidget(self.live_results_table)
            
            # Export button
            export_layout = QHBoxLayout()
            export_csv_btn = QPushButton("📄 Export CSV")
            export_csv_btn.clicked.connect(self.export_live_results_csv)
            export_csv_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; padding: 8px;")
            export_layout.addWidget(export_csv_btn)
            export_layout.addStretch()
            tab_layout.addLayout(export_layout)
            
            self.live_results_tab.setLayout(tab_layout)
            self.tab_widget_ref.addTab(self.live_results_tab, tab_name)
            tab_index = self.tab_widget_ref.count() - 1
        
        # Update info label
        num_isotopes = len(self.live_batch_results)
        total_lines = sum(len(data.get('energies', [])) for data in self.live_batch_results.values())
        self.live_info.setText(
            f"<b>Completed Isotopes:</b> {num_isotopes} | "
            f"<b>Energy Lines:</b> {total_lines} | "
            f"<b>BeamOn:</b> {beamon:,}"
        )
        
        # Clear and repopulate table
        self.live_results_table.setRowCount(0)
        
        row_idx = 0
        for isotope, data in self.live_batch_results.items():
            peaks = ISOTOPE_DATABASE.get(isotope, [])
            num_energies = len(peaks)
            has_cascade = isotope in cascade_isotopes
            
            energy_results = data.get('energies', [])
            
            for e_data in energy_results:
                self.live_results_table.insertRow(row_idx)
                
                energy = e_data['energy']
                intensity = e_data['intensity']
                cnt_a = e_data['cnt_a']
                cnt_b = e_data['cnt_b']
                
                # Calculate emissions
                if num_energies == 1:
                    emissions = beamon
                else:
                    emissions = beamon * (intensity / 100.0)
                
                # Mode B correction: Custom mode divides beamon among energies
                # So cnt_b must be multiplied by num_energies for correct efficiency
                if num_energies > 1:
                    cnt_b_adj = cnt_b * num_energies
                else:
                    cnt_b_adj = cnt_b
                
                # Efficiency calculations
                eff_a = cnt_a / emissions if emissions > 0 else 0
                eff_b = cnt_b_adj / emissions if emissions > 0 else 0  # Use adjusted count!
                
                # COI calculation
                coi = cnt_a / cnt_b_adj if cnt_b_adj > 0 else 1.0
                sigma_coi = coi * np.sqrt(1.0/cnt_a + 1.0/cnt_b_adj) if cnt_a > 0 and cnt_b_adj > 0 else 0
                
                # Populate row
                isotope_item = QTableWidgetItem(isotope)
                if has_cascade:
                    isotope_item.setToolTip("Cascade isotope")
                    isotope_item.setForeground(QColor("#E65100"))
                self.live_results_table.setItem(row_idx, 0, isotope_item)
                
                self.live_results_table.setItem(row_idx, 1, QTableWidgetItem(f"{energy:.3f}"))
                self.live_results_table.setItem(row_idx, 2, QTableWidgetItem(f"{intensity:.2f}"))
                self.live_results_table.setItem(row_idx, 3, QTableWidgetItem(f"{cnt_a:,}"))
                self.live_results_table.setItem(row_idx, 4, QTableWidgetItem(f"{cnt_b:,}"))
                self.live_results_table.setItem(row_idx, 5, QTableWidgetItem(f"{eff_a:.6f}"))
                self.live_results_table.setItem(row_idx, 6, QTableWidgetItem(f"{eff_b:.6f}"))
                
                # COI with color coding
                coi_item = QTableWidgetItem(f"{coi:.4f}")
                if has_cascade and abs(coi - 1.0) > 0.01:
                    coi_item.setBackground(QColor("#fff3cd"))  # Yellow for cascade effect
                elif not has_cascade and abs(coi - 1.0) > 0.05:
                    coi_item.setBackground(QColor("#f8d7da"))  # Red for unexpected
                else:
                    coi_item.setBackground(QColor("#d4edda"))  # Green for OK
                self.live_results_table.setItem(row_idx, 7, coi_item)
                
                self.live_results_table.setItem(row_idx, 8, QTableWidgetItem(f"{coi:.4f} ± {sigma_coi:.4f}"))
                
                row_idx += 1
        
        # Switch to results tab
        self.tab_widget_ref.setCurrentIndex(tab_index)

    def export_live_results_csv(self):
        """Export live results to CSV"""
        if not self.live_batch_results:
            QMessageBox.warning(self, "No Data", "No results to export.")
            return
        
        filename, _ = QFileDialog.getSaveFileName(self, "Export Live Results", 
                                                   "live_batch_results.csv", "CSV (*.csv)")
        if not filename:
            return
        
        beamon = self.sim_widget_ref.num_events.value() if self.sim_widget_ref else 100000
        cascade_isotopes = ['Co60', 'Y88', 'Na22', 'Eu152', 'Ba133']
        
        try:
            with open(filename, 'w') as f:
                f.write("Isotope,Energy_keV,Intensity_pct,Counts_A,Counts_B,Eff_A,Eff_B,COI,Sigma_COI,Has_Cascade\n")
                
                for isotope, data in self.live_batch_results.items():
                    peaks = ISOTOPE_DATABASE.get(isotope, [])
                    num_energies = len(peaks)
                    has_cascade = isotope in cascade_isotopes
                    
                    for e_data in data.get('energies', []):
                        energy = e_data['energy']
                        intensity = e_data['intensity']
                        cnt_a = e_data['cnt_a']
                        cnt_b = e_data['cnt_b']
                        
                        if num_energies == 1:
                            emissions = beamon
                        else:
                            emissions = beamon * (intensity / 100.0)
                        
                        # Mode B correction: Custom mode divides beamon among energies
                        if num_energies > 1:
                            cnt_b_adj = cnt_b * num_energies
                        else:
                            cnt_b_adj = cnt_b
                        
                        eff_a = cnt_a / emissions if emissions > 0 else 0
                        eff_b = cnt_b_adj / emissions if emissions > 0 else 0  # Use adjusted count!
                        
                        coi = cnt_a / cnt_b_adj if cnt_b_adj > 0 else 1.0
                        sigma_coi = coi * np.sqrt(1.0/cnt_a + 1.0/cnt_b_adj) if cnt_a > 0 and cnt_b_adj > 0 else 0
                        
                        f.write(f"{isotope},{energy:.3f},{intensity:.2f},{cnt_a},{cnt_b},"
                               f"{eff_a:.8f},{eff_b:.8f},{coi:.6f},{sigma_coi:.6f},{has_cascade}\n")
            
            QMessageBox.information(self, "Export", f"Results exported to:\n{filename}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Error exporting:\n{str(e)}")

    def clear_live_results(self):
        """Clear live results (called when starting a new batch)"""
        self.live_batch_results = {}

    def load_results_dialog(self):
        default_path = os.path.join(self.sim_widget_ref.sim_path_edit.text(), "output.dat") if self.sim_widget_ref else ""
        filename, _ = QFileDialog.getOpenFileName(self, "Open Result", default_path, "Data (*.dat *.txt)")
        if filename: self.load_results_file(filename)

    def load_results_file(self, filename):
        try:
            data = np.loadtxt(filename, comments='#')
            if len(data.shape) < 2: return
            
            energies = []
            for row in data:
                if len(row) < 3: continue
                n_det = int(row[2])
                for i in range(n_det):
                    idx = 4 + i*3
                    if idx < len(row): energies.append(row[idx])
            
            self.energies = np.array(energies)
            self.canvas.plot_spectrum(self.energies, self.current_isotope)
            self.populate_peaks_table()
            self.stats_text.setText(f"Loaded {len(energies)} energy deposits from {filename}.")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load: {e}")

    def on_isotope_changed(self, txt):
        self.current_isotope = txt
        if self.energies is not None: self.canvas.plot_spectrum(self.energies, txt)

    def populate_peaks_table(self):
        self.peaks_table.setRowCount(0)
        if not hasattr(self.canvas, 'detected_peaks'): return
        
        for p in self.canvas.detected_peaks:
            r = self.peaks_table.rowCount()
            self.peaks_table.insertRow(r)
            self.peaks_table.setItem(r,0, QTableWidgetItem(p['name']))
            self.peaks_table.setItem(r,1, QTableWidgetItem(f"{p['theoretical_energy']:.2f}"))
            self.peaks_table.setItem(r,2, QTableWidgetItem(f"{p['energy']:.2f}"))
            
            w_l, w_h = p['energy']*0.97, p['energy']*1.03
            cnts = np.sum((self.energies >= w_l) & (self.energies <= w_h))
            self.peaks_table.setItem(r,3, QTableWidgetItem(f"{int(cnts)}"))
            self.peaks_table.setItem(r,4, QTableWidgetItem(f"{int(cnts)}"))
            self.peaks_table.setItem(r,5, QTableWidgetItem(f"{p['intensity']:.2f}"))
            
            emit = self.beamon_input.value() * (p['intensity']/100)
            eff = cnts/emit if emit > 0 else 0
            self.peaks_table.setItem(r,6, QTableWidgetItem(f"{eff:.4e}"))

    def on_peak_selected(self, e, c):
        QMessageBox.information(self, "Peak", f"Energy: {e:.2f} keV\nCounts: {c}")

    def calculate_and_show_coi(self):
        if not self.sim_widget_ref:
             QMessageBox.critical(self, "Error", "Simulation Control reference is missing.")
             return
             
        sim_path = self.sim_widget_ref.sim_path_edit.text()
        file_a = os.path.join(sim_path, "output_mode_a_tcs.txt")
        file_b = os.path.join(sim_path, "output_mode_b_notcs.txt")
        
        if not (os.path.exists(file_a) and os.path.exists(file_b)):
            QMessageBox.warning(self, "Missing Files", 
                "Could not find Mode A/B output files.\nPlease run the batch analysis first.")
            return

        try:
            energies_a = self.load_raw_energies(file_a)
            energies_b = self.load_raw_energies(file_b)
            
            self.show_coi_results_tab(energies_a, energies_b)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Analysis failed: {e}")

    def load_multi_isotope_files(self):
        """Load multi-isotope A/B files from batch run and display results"""
        if not self.sim_widget_ref:
            QMessageBox.critical(self, "Error", "Simulation Control reference is missing.")
            return
        
        sim_path = self.sim_widget_ref.sim_path_edit.text()
        
        multi_isotope_results = {}
        found_files = []
        
        for isotope in ISOTOPE_DATABASE.keys():
            file_a = os.path.join(sim_path, f"output_{isotope}_mode_a.txt")
            file_b = os.path.join(sim_path, f"output_{isotope}_mode_b.txt")
            
            if os.path.exists(file_a) and os.path.exists(file_b):
                found_files.append(isotope)
                try:
                    energies_a = self.load_raw_energies(file_a)
                    energies_b = self.load_raw_energies(file_b)
                    
                    isotope_peaks = ISOTOPE_DATABASE.get(isotope, [])
                    energy_results = []
                    
                    for peak in isotope_peaks:
                        energy = peak['energy']
                        intensity = peak['intensity']
                        
                        cnt_a = self.get_peak_counts(energies_a, energy)
                        cnt_b = self.get_peak_counts(energies_b, energy)
                        
                        energy_results.append({
                            'energy': energy,
                            'intensity': intensity,
                            'cnt_a': int(cnt_a),
                            'cnt_b': int(cnt_b)
                        })
                    
                    multi_isotope_results[isotope] = {'energies': energy_results}
                    
                except Exception as e:
                    print(f"Error loading {isotope}: {e}")
        
        if not found_files:
            QMessageBox.warning(self, "No Files Found", 
                "No multi-isotope output files found.\n\n"
                "Expected files: output_ISOTOPE_mode_a.txt and output_ISOTOPE_mode_b.txt\n"
                "Run 'Multi-Isotope Batch' first to generate these files.")
            return
        
        msg = f"Found {len(found_files)} isotope(s) with A/B files:\n\n"
        msg += ", ".join(found_files)
        msg += "\n\nLoad and analyze these results?"
        
        reply = QMessageBox.question(self, "Multi-Isotope Files Found", msg,
                                    QMessageBox.Yes | QMessageBox.No)
        
        if reply == QMessageBox.Yes:
            beamon = self.sim_widget_ref.num_events.value()
            
            self.batch_dialog = BatchResultsDialog(self, multi_isotope_results, beamon)
            self.batch_dialog.show()

    def load_raw_energies(self, path):
        data = np.loadtxt(path, comments='#')
        energies = []
        if len(data.shape) < 2: return np.array(energies)
        
        for row in data:
            if len(row) < 3: continue
            n_det = int(row[2])
            for i in range(n_det):
                idx = 4 + i*3
                if idx < len(row): energies.append(row[idx])
        return np.array(energies)

    def get_peak_counts(self, energies, peak_energy):
        """Get counts in 3% window around peak energy"""
        w_l = peak_energy * 0.97
        w_h = peak_energy * 1.03
        counts = np.sum((energies >= w_l) & (energies <= w_h))
        return counts

    def show_coi_results_tab(self, energies_a, energies_b):
        tab_name = "COI Results"
        idx = -1
        for i in range(self.tab_widget_ref.count()):
            if self.tab_widget_ref.tabText(i) == tab_name: idx = i; break
        
        content = QWidget()
        layout = QVBoxLayout()
        
        iso = self.config_widget_ref.single_isotope_combo.currentText()
        isotope_peaks = ISOTOPE_DATABASE.get(iso, [])
        num_energies = len(isotope_peaks)
        beamon = self.beamon_input.value()
        
        cascade_isotopes = ['Co60', 'Y88', 'Na22', 'Eu152', 'Ba133']
        has_cascade = iso in cascade_isotopes
        
        # Create results table
        table = QTableWidget()
        table.setColumnCount(7)
        table.setHorizontalHeaderLabels([
            "Energy (keV)", "Counts A", "Counts B", "Eff A", "Eff B", "COI", "COI ± σ"
        ])
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        for i, peak in enumerate(isotope_peaks):
            t_energy = peak['energy']
            intensity = peak['intensity']
            
            cnt_a = self.get_peak_counts(energies_a, t_energy)
            cnt_b_raw = self.get_peak_counts(energies_b, t_energy)
            
            # Mode B correction: Custom mode divides beamon among energies
            # Apply for ALL multi-energy isotopes, not just cascade
            if num_energies > 1:
                cnt_b = num_energies * cnt_b_raw
            else:
                cnt_b = cnt_b_raw
            
            if num_energies == 1:
                emissions = beamon
            else:
                emissions = beamon * (intensity / 100.0)
            
            eff_a = (cnt_a / emissions) if emissions > 0 else 0.0
            eff_b = (cnt_b / emissions) if emissions > 0 else 0.0
            
            if has_cascade:
                coi = cnt_a / cnt_b if cnt_b > 0 else 1.0
            else:
                coi = cnt_a / cnt_b if cnt_b > 0 else 1.0
            
            if cnt_a > 0 and cnt_b > 0:
                sigma_coi = coi * np.sqrt(1.0/cnt_a + 1.0/cnt_b)
            else:
                sigma_coi = 0.0
            
            table.insertRow(i)
            table.setItem(i, 0, QTableWidgetItem(f"{t_energy:.3f}"))
            table.setItem(i, 1, QTableWidgetItem(str(int(cnt_a))))
            table.setItem(i, 2, QTableWidgetItem(str(int(cnt_b))))
            table.setItem(i, 3, QTableWidgetItem(f"{eff_a:.6f}"))
            table.setItem(i, 4, QTableWidgetItem(f"{eff_b:.6f}"))
            
            coi_item = QTableWidgetItem(f"{coi:.4f}")
            if has_cascade and abs(coi - 1.0) > 0.01:
                coi_item.setBackground(QColor("#fff3cd"))
            elif not has_cascade and abs(coi - 1.0) > 0.05:
                coi_item.setBackground(QColor("#f8d7da"))
            table.setItem(i, 5, coi_item)
            table.setItem(i, 6, QTableWidgetItem(f"{coi:.4f} ± {sigma_coi:.4f}"))
        
        layout.addWidget(table)
        content.setLayout(layout)
        
        if idx != -1: self.tab_widget_ref.removeTab(idx)
        self.tab_widget_ref.addTab(content, tab_name)
        self.tab_widget_ref.setCurrentIndex(self.tab_widget_ref.count()-1)



# ==================================================================
# GEOMETRY OPTIMIZATION WIDGET v3.0
# Follows EXACTLY the same logic as "Run Multi-Isotope Batch Analysis"
# - Mode A: Full isotope with TCS (cascade emission)
# - Mode B: Custom isotope without TCS (independent emission)
# - Efficiency calculation: counts / (beamon × intensity/100)
# ==================================================================

class OptimizationWorker(QThread):
    """
    Worker for geometry optimization.
    Follows SAME logic as Run Multi-Isotope Batch Analysis.
    """
    progress = pyqtSignal(int, int, str)
    iteration_result = pyqtSignal(dict)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    log = pyqtSignal(str)
    
    def __init__(self, config):
        super().__init__()
        self.config = config
        self._is_running = True
        
    def stop(self):
        self._is_running = False
        
    def run(self):
        try:
            param_name = self.config['param_name']
            start = self.config['start']
            stop = self.config['stop']
            step = self.config['step']
            mode = self.config.get('mode', 'A')
            
            values = np.arange(start, stop + step/2, step)
            n_values = len(values)
            
            # Get isotopes and their energies
            isotopes_data = self.config['isotopes_data']  # {isotope: [{energy, intensity, exp_eff, exp_unc}]}
            
            results = []
            best_chi2 = float('inf')
            best_value = start
            best_result = None
            
            self.log.emit(f"\n{'='*60}")
            self.log.emit(f"GEOMETRY OPTIMIZATION v3.0")
            self.log.emit(f"Following Multi-Isotope Batch Logic")
            self.log.emit(f"{'='*60}")
            self.log.emit(f"Parameter: {param_name}")
            self.log.emit(f"Range: {start:.3f} to {stop:.3f} mm, step {step:.3f}")
            self.log.emit(f"Mode: {'A (TCS/cascade)' if mode == 'A' else 'B (independent)'}")
            self.log.emit(f"Parameter iterations: {n_values}")
            self.log.emit(f"Isotopes: {list(isotopes_data.keys())}")
            self.log.emit(f"Events (beamon): {self.config['num_events']}")
            self.log.emit(f"{'='*60}\n")
            
            for i_val, value in enumerate(values):
                if not self._is_running:
                    self.log.emit("\n⚠️ Optimization stopped by user")
                    break
                
                self.progress.emit(i_val + 1, n_values, 
                                  f"{param_name} = {value:.4f} mm")
                
                self.log.emit(f"\n[{i_val+1}/{n_values}] {param_name} = {value:.4f} mm")
                self.log.emit("-" * 50)
                
                # Results for this configuration
                all_efficiencies = {}
                all_counts = {}
                
                # Process each isotope
                for isotope, peaks in isotopes_data.items():
                    self.log.emit(f"\n  Isotope: {isotope}")
                    
                    energies = [p['energy'] for p in peaks]
                    intensities = [p['intensity'] for p in peaks]
                    
                    # Generate config for this isotope
                    self.generate_isotope_config(param_name, value, isotope, 
                                                energies, intensities, mode)
                    
                    # Run simulation
                    output_path = self.run_simulation()
                    
                    if output_path is None:
                        self.log.emit(f"    ⚠️ Simulation failed")
                        continue
                    
                    # Load raw energies from output
                    raw_energies = self.load_raw_energies(output_path)
                    self.log.emit(f"    Total detected events: {len(raw_energies)}")
                    
                    # Count peaks and calculate efficiency (same as Batch)
                    num_energies = len(peaks)
                    beamon = self.config['num_events']
                    
                    for peak in peaks:
                        energy = peak['energy']
                        intensity = peak['intensity']
                        exp_eff = peak.get('exp_eff', 0)
                        exp_unc = peak.get('exp_unc', exp_eff * 0.05)
                        
                        # Count in ±3% window
                        counts_raw = self.get_peak_counts(raw_energies, energy)
                        
                        # Mode B correction: Custom mode divides beamon among energies
                        # So counts must be multiplied by num_energies for correct efficiency
                        if mode == 'B' and num_energies > 1:
                            counts = counts_raw * num_energies
                        else:
                            counts = counts_raw
                        
                        # Calculate emissions (SAME as Batch)
                        if num_energies == 1:
                            emissions = beamon
                        else:
                            emissions = beamon * (intensity / 100.0)
                        
                        # Calculate efficiency
                        eff = counts / emissions if emissions > 0 else 0
                        
                        all_efficiencies[energy] = eff
                        all_counts[energy] = counts
                        
                        ratio = eff / exp_eff if exp_eff > 0 else 0
                        if mode == 'B' and num_energies > 1:
                            self.log.emit(f"    {energy:.1f} keV: cnt_raw={counts_raw}, cnt_adj={counts}, "
                                         f"eff={eff:.6f}, exp={exp_eff:.6f}, ratio={ratio:.3f}")
                        else:
                            self.log.emit(f"    {energy:.1f} keV: cnt={counts}, "
                                         f"eff={eff:.6f}, exp={exp_eff:.6f}, ratio={ratio:.3f}")
                
                # Calculate chi2
                chi2 = self.calculate_chi2(all_efficiencies, isotopes_data)
                
                result = {
                    'param_value': value,
                    'chi2': chi2,
                    'efficiencies': all_efficiencies,
                    'counts': all_counts
                }
                results.append(result)
                
                self.log.emit(f"\n  → χ² reduced = {chi2:.4f}")
                
                if chi2 < best_chi2:
                    best_chi2 = chi2
                    best_value = value
                    best_result = result
                    self.log.emit(f"  ★ New best!")
                
                self.iteration_result.emit(result)
            
            final_results = {
                'param_name': param_name,
                'all_results': results,
                'best_value': best_value,
                'best_chi2': best_chi2,
                'best_result': best_result,
                'mode': mode
            }
            
            self.log.emit(f"\n{'='*60}")
            self.log.emit(f"OPTIMIZATION COMPLETE")
            self.log.emit(f"Best {param_name} = {best_value:.4f} mm")
            self.log.emit(f"Best χ² = {best_chi2:.4f}")
            self.log.emit(f"{'='*60}\n")
            
            self.finished.emit(final_results)
            
        except Exception as e:
            import traceback
            self.log.emit(f"ERROR: {str(e)}")
            self.log.emit(traceback.format_exc())
            self.error.emit(str(e))
    
    def generate_isotope_config(self, param_name, param_value, isotope, 
                                energies, intensities, mode):
        """
        Generate config.txt following Multi-Isotope Batch logic:
        - Mode A: Use isotope name directly, enable_tcs=true
        - Mode B: Use Custom with energies/intensities, enable_tcs=false
        """
        config_path = os.path.join(self.config['sim_path'], 'config.txt')
        det_params = self.config['detector_params']
        source_config = self.config.get('source_config', {})
        
        # Parameter mapping
        param_mapping = {
            'crystal_diameter': 'crystal_diameter',
            'crystal_length': 'crystal_length',
            'hole_diameter': 'hole_diameter',
            'hole_depth': 'hole_depth',
            'front_dead_layer': 'ge_dead_layer_front',
            'lateral_dead_layer': 'ge_dead_layer',
            'li_dead_layer': 'li_dead_layer',
            'li_dead_layer_front': 'li_dead_layer_front',
            'al_window_thickness': 'al_window_thickness',
            'al_cup_thickness': 'al_cup_thickness',
            'vacuum_gap': 'vacuum_gap'
        }
        
        lines = []
        
        # Detector geometry
        lines.append("# === DETECTOR GEOMETRY ===")
        for gui_name, config_name in param_mapping.items():
            if gui_name == param_name:
                lines.append(f"{config_name} = {param_value}")
            else:
                lines.append(f"{config_name} = {det_params.get(gui_name, 0)}")
        
        lines.append("")
        lines.append("# === ISOTOPE CONFIGURATION ===")
        
        if mode == 'A':
            # Mode A: Use isotope name directly (cascade emission)
            lines.append(f"isotope = {isotope}")
            lines.append("# Mode A: Full isotope with cascade emission")
        else:
            # Mode B: Use Custom (independent emission)
            energies_str = ",".join([f"{e:.2f}" for e in energies])
            intensities_str = ",".join([f"{i:.2f}" for i in intensities])
            lines.append("isotope = Custom")
            lines.append(f"custom_energies = {energies_str}")
            lines.append(f"custom_intensities = {intensities_str}")
            lines.append("# Mode B: Independent single-gamma emission")
        
        # Source geometry
        lines.append("")
        lines.append("# === SOURCE GEOMETRY ===")
        source_type = source_config.get('source_type', 'point')
        lines.append(f"source_type = {source_type}")
        
        if source_type == 'marinelli':
            lines.append(f"marinelli_type = {source_config.get('marinelli_type', source_config.get('type', '1000ml'))}")
            lines.append(f"marinelli_fill_material = {source_config.get('marinelli_fill_material', source_config.get('fill_material', 'water'))}")
        elif source_type == 'disk':
            lines.append(f"disk_radius = {source_config.get('disk_radius', source_config.get('radius', 25.0))}")
            lines.append(f"disk_thickness = {source_config.get('disk_thickness', source_config.get('thickness', 5.0))}")
            lines.append(f"disk_position_z = {source_config.get('disk_position_z', source_config.get('position_z', -50.0))}")
            lines.append(f"disk_material = {source_config.get('disk_material', source_config.get('material', 'water'))}")
        elif source_type == 'volume':
            lines.append(f"cylinder_inner_radius = {source_config.get('cylinder_inner_radius', source_config.get('inner_radius', 28.0))}")
            lines.append(f"cylinder_outer_radius = {source_config.get('cylinder_outer_radius', source_config.get('outer_radius', 30.0))}")
            lines.append(f"cylinder_height = {source_config.get('cylinder_height', source_config.get('height', 50.0))}")
            lines.append(f"cylinder_wall_thickness = {source_config.get('cylinder_wall_thickness', source_config.get('wall_thickness', 2.0))}")
            lines.append(f"cylinder_bottom_thickness = {source_config.get('cylinder_bottom_thickness', source_config.get('bottom_thickness', 2.0))}")
            lines.append(f"cylinder_position_z = {source_config.get('cylinder_position_z', source_config.get('position_z', -60.0))}")
            lines.append(f"cylinder_wall_material = {source_config.get('cylinder_wall_material', source_config.get('wall_material', 'polypropylene'))}")
            lines.append(f"cylinder_fill_material = {source_config.get('cylinder_fill_material', source_config.get('fill_material', 'water'))}")
        else:
            lines.append(f"source_x = {source_config.get('source_x', 0.0)}")
            lines.append(f"source_y = {source_config.get('source_y', 0.0)}")
            lines.append(f"source_z = {source_config.get('source_z', 5.0)}")
        
        # Simulation settings
        lines.append("")
        lines.append("# === SIMULATION SETTINGS ===")
        lines.append("energy_resolution = 2.0")
        lines.append("energy_threshold = 10.0")
        lines.append("coincidence_window = 500.0")
        
        # TCS setting based on mode
        if mode == 'A':
            lines.append("enable_tcs = true")
        else:
            lines.append("enable_tcs = false")
        
        lines.append("fast_mode = true")
        
        with open(config_path, 'w') as f:
            f.write('\n'.join(lines))
    
    def run_simulation(self):
        """Run simulation and return output path"""
        sim_path = self.config['sim_path']
        geant4_path = self.config['geant4_path']
        num_events = self.config['num_events']
        
        # Create macro
        macro_path = os.path.join(sim_path, 'run_opt.mac')
        with open(macro_path, 'w') as f:
            f.write(f'/run/initialize\n/run/beamOn {num_events}\n')
        
        # Create script
        script_path = os.path.join(sim_path, 'run_opt.sh')
        with open(script_path, 'w') as f:
            f.write(f'#!/bin/bash\n')
            f.write(f'source {geant4_path}\n')
            f.write(f'cd {sim_path}\n')
            f.write(f'./hpge_sim run_opt.mac 2>&1\n')
        os.chmod(script_path, 0o755)
        
        try:
            process = subprocess.Popen(
                ['bash', script_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                cwd=sim_path
            )
            stdout, _ = process.communicate(timeout=1800)
            
            if process.returncode != 0:
                return None
            
            output_path = os.path.join(sim_path, 'output.dat')
            if os.path.exists(output_path):
                return output_path
            return None
            
        except subprocess.TimeoutExpired:
            process.kill()
            return None
        except Exception as e:
            self.log.emit(f"    Error: {str(e)}")
            return None
    
    def load_raw_energies(self, path):
        """Load raw energies from output file (SAME as Batch)"""
        data = np.loadtxt(path, comments='#')
        energies = []
        if len(data.shape) < 2:
            return np.array(energies)
        
        for row in data:
            if len(row) < 3:
                continue
            n_det = int(row[2])
            for i in range(n_det):
                idx = 4 + i * 3
                if idx < len(row):
                    energies.append(row[idx])
        
        return np.array(energies)
    
    def get_peak_counts(self, energies, peak_energy):
        """Get counts in ±3% window around peak (SAME as Batch)"""
        w_l = peak_energy * 0.97
        w_h = peak_energy * 1.03
        counts = np.sum((energies >= w_l) & (energies <= w_h))
        return int(counts)
    
    def calculate_chi2(self, sim_efficiencies, isotopes_data):
        """Calculate reduced chi-squared"""
        chi2 = 0.0
        n_points = 0
        
        for isotope, peaks in isotopes_data.items():
            for peak in peaks:
                energy = peak['energy']
                exp_eff = peak.get('exp_eff', 0)
                exp_unc = peak.get('exp_unc', exp_eff * 0.05)
                
                if energy in sim_efficiencies and exp_unc > 0 and exp_eff > 0:
                    sim_eff = sim_efficiencies[energy]
                    chi2 += ((sim_eff - exp_eff) / exp_unc) ** 2
                    n_points += 1
        
        if n_points > 1:
            return chi2 / (n_points - 1)
        return chi2


class IsotopeEfficiencyDialog(QDialog):
    """Dialog to select isotopes and enter experimental efficiencies"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Configure Isotopes & Experimental Efficiencies")
        self.setMinimumSize(900, 700)
        self.isotopes_data = {}
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        
        # Instructions
        info = QLabel(
            "<b>Instructions:</b><br>"
            "1. Select isotopes from the left panel<br>"
            "2. Enter experimental efficiencies in the right table<br>"
            "3. Efficiencies are required for χ² calculation"
        )
        info.setStyleSheet("background: #E3F2FD; padding: 10px; border-radius: 5px;")
        layout.addWidget(info)
        
        splitter = QSplitter(Qt.Horizontal)
        
        # === LEFT: Isotope selection ===
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        
        left_layout.addWidget(QLabel("<b>Select Isotopes:</b>"))
        
        self.isotope_list = QListWidget()
        self.isotope_list.setSelectionMode(QListWidget.MultiSelection)
        
        cascade_isotopes = ['Co60', 'Y88', 'Na22', 'Eu152', 'Ba133']
        
        for iso_name, peaks in ISOTOPE_DATABASE.items():
            energies_str = ", ".join([f"{p['energy']:.1f}" for p in peaks])
            cascade = " ⚡CASCADE" if iso_name in cascade_isotopes else ""
            item = QListWidgetItem(f"{iso_name}{cascade}\n  {energies_str} keV")
            item.setData(Qt.UserRole, iso_name)
            if iso_name in cascade_isotopes:
                item.setForeground(QColor("#D84315"))
            self.isotope_list.addItem(item)
        
        left_layout.addWidget(self.isotope_list)
        
        btn_layout = QHBoxLayout()
        select_all = QPushButton("Select All")
        select_all.clicked.connect(self.isotope_list.selectAll)
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.isotope_list.clearSelection)
        btn_layout.addWidget(select_all)
        btn_layout.addWidget(clear_btn)
        left_layout.addLayout(btn_layout)
        
        load_btn = QPushButton("→ Load Selected to Table")
        load_btn.setStyleSheet("background: #4CAF50; color: white; font-weight: bold; padding: 10px;")
        load_btn.clicked.connect(self.load_to_table)
        left_layout.addWidget(load_btn)
        
        splitter.addWidget(left_widget)
        
        # === RIGHT: Efficiency table ===
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        
        right_layout.addWidget(QLabel("<b>Experimental Efficiencies:</b>"))
        
        csv_layout = QHBoxLayout()
        csv_btn = QPushButton("📁 Load from CSV")
        csv_btn.clicked.connect(self.load_csv)
        csv_layout.addWidget(csv_btn)
        csv_layout.addStretch()
        right_layout.addLayout(csv_layout)
        
        self.eff_table = QTableWidget()
        self.eff_table.setColumnCount(5)
        self.eff_table.setHorizontalHeaderLabels([
            "Isotope", "Energy (keV)", "Intensity (%)", "Exp. Efficiency", "Uncertainty"
        ])
        self.eff_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        right_layout.addWidget(self.eff_table)
        
        splitter.addWidget(right_widget)
        splitter.setSizes([300, 600])
        
        layout.addWidget(splitter)
        
        # Dialog buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.validate_and_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
        self.setLayout(layout)
    
    def load_to_table(self):
        """Load selected isotopes to the efficiency table"""
        selected = []
        for item in self.isotope_list.selectedItems():
            selected.append(item.data(Qt.UserRole))
        
        if not selected:
            QMessageBox.warning(self, "Warning", "Select at least one isotope")
            return
        
        # Preserve existing data
        existing = {}
        for row in range(self.eff_table.rowCount()):
            try:
                e_item = self.eff_table.item(row, 1)
                eff_item = self.eff_table.item(row, 3)
                unc_item = self.eff_table.item(row, 4)
                if e_item and e_item.text():
                    e = float(e_item.text())
                    existing[e] = {
                        'eff': float(eff_item.text()) if eff_item and eff_item.text() else 0,
                        'unc': float(unc_item.text()) if unc_item and unc_item.text() else 0
                    }
            except:
                continue
        
        # Build table
        rows = []
        for iso in selected:
            for peak in ISOTOPE_DATABASE.get(iso, []):
                rows.append({
                    'isotope': iso,
                    'energy': peak['energy'],
                    'intensity': peak['intensity']
                })
        
        rows.sort(key=lambda x: x['energy'])
        
        self.eff_table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            # Isotope (read-only)
            iso_item = QTableWidgetItem(row['isotope'])
            iso_item.setFlags(iso_item.flags() & ~Qt.ItemIsEditable)
            self.eff_table.setItem(i, 0, iso_item)
            
            # Energy (read-only)
            e_item = QTableWidgetItem(f"{row['energy']:.2f}")
            e_item.setFlags(e_item.flags() & ~Qt.ItemIsEditable)
            self.eff_table.setItem(i, 1, e_item)
            
            # Intensity (read-only)
            int_item = QTableWidgetItem(f"{row['intensity']:.1f}")
            int_item.setFlags(int_item.flags() & ~Qt.ItemIsEditable)
            self.eff_table.setItem(i, 2, int_item)
            
            # Efficiency (editable) - try to find existing
            eff_val = ""
            unc_val = ""
            for exist_e, exist_data in existing.items():
                if abs(exist_e - row['energy']) < 1:
                    eff_val = f"{exist_data['eff']:.6g}" if exist_data['eff'] > 0 else ""
                    unc_val = f"{exist_data['unc']:.6g}" if exist_data['unc'] > 0 else ""
                    break
            
            self.eff_table.setItem(i, 3, QTableWidgetItem(eff_val))
            self.eff_table.setItem(i, 4, QTableWidgetItem(unc_val))
        
        QMessageBox.information(self, "Loaded", 
                               f"Loaded {len(rows)} energy lines from {len(selected)} isotopes.\n"
                               f"Please enter experimental efficiencies.")
    
    def load_csv(self):
        """Load experimental efficiencies from CSV"""
        filename, _ = QFileDialog.getOpenFileName(self, "Load CSV", "", 
                                                   "CSV (*.csv);;All (*)")
        if not filename:
            return
        
        try:
            csv_data = []
            with open(filename, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.replace(';', ',').split(',')
                    if len(parts) >= 2:
                        try:
                            e = float(parts[0])
                            eff = float(parts[1])
                            unc = float(parts[2]) if len(parts) >= 3 else eff * 0.05
                            csv_data.append((e, eff, unc))
                        except:
                            continue
            
            # Match with table
            matched = 0
            for csv_e, csv_eff, csv_unc in csv_data:
                for row in range(self.eff_table.rowCount()):
                    e_item = self.eff_table.item(row, 1)
                    if e_item:
                        table_e = float(e_item.text())
                        if abs(table_e - csv_e) < 2:  # Within 2 keV
                            self.eff_table.setItem(row, 3, QTableWidgetItem(f"{csv_eff:.6g}"))
                            self.eff_table.setItem(row, 4, QTableWidgetItem(f"{csv_unc:.6g}"))
                            matched += 1
                            break
            
            QMessageBox.information(self, "Loaded", 
                                   f"Matched {matched} efficiencies from CSV")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
    
    def validate_and_accept(self):
        """Validate data and accept"""
        self.isotopes_data = {}
        
        for row in range(self.eff_table.rowCount()):
            try:
                iso = self.eff_table.item(row, 0).text()
                energy = float(self.eff_table.item(row, 1).text())
                intensity = float(self.eff_table.item(row, 2).text())
                
                eff_item = self.eff_table.item(row, 3)
                unc_item = self.eff_table.item(row, 4)
                
                exp_eff = float(eff_item.text()) if eff_item and eff_item.text() else 0
                exp_unc = float(unc_item.text()) if unc_item and unc_item.text() else exp_eff * 0.05
                
                if exp_eff > 0:  # Only include if efficiency provided
                    if iso not in self.isotopes_data:
                        self.isotopes_data[iso] = []
                    
                    self.isotopes_data[iso].append({
                        'energy': energy,
                        'intensity': intensity,
                        'exp_eff': exp_eff,
                        'exp_unc': exp_unc
                    })
            except:
                continue
        
        if not self.isotopes_data:
            QMessageBox.warning(self, "Error", 
                               "No valid experimental efficiencies entered!\n"
                               "Please enter at least one efficiency value.")
            return
        
        total_energies = sum(len(v) for v in self.isotopes_data.values())
        self.accept()
    
    def get_data(self):
        return self.isotopes_data


class GeometryOptimizationWidget(QWidget):
    """
    Geometry Optimization Widget v3.0
    Follows EXACTLY the same logic as Run Multi-Isotope Batch Analysis
    """
    
    def __init__(self):
        super().__init__()
        self.isotopes_data = {}
        self.optimization_results = []
        self.worker = None
        self.sim_widget_ref = None
        self.config_widget_ref = None
        self.best_value = None
        self.best_param = None
        self.init_ui()
    
    def set_references(self, config_widget, sim_widget):
        self.config_widget_ref = config_widget
        self.sim_widget_ref = sim_widget
    
    def init_ui(self):
        main_layout = QVBoxLayout()
        
        # Title
        title = QLabel("🔧 Detector Geometry Optimization v3.0")
        title.setFont(QFont("Arial", 14, QFont.Bold))
        title.setStyleSheet("color: #1565C0;")
        main_layout.addWidget(title)
        
        subtitle = QLabel("Same logic as 'Run Multi-Isotope Batch Analysis'")
        subtitle.setStyleSheet("color: #666; font-style: italic;")
        main_layout.addWidget(subtitle)
        
        splitter = QSplitter(Qt.Horizontal)
        
        # === LEFT PANEL ===
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # 1. Isotope Configuration
        iso_group = QGroupBox("1. Isotopes & Experimental Efficiencies")
        iso_layout = QVBoxLayout()
        
        self.config_iso_btn = QPushButton("📊 Configure Isotopes && Efficiencies")
        self.config_iso_btn.setStyleSheet(
            "background: #2196F3; color: white; font-weight: bold; padding: 10px;"
        )
        self.config_iso_btn.clicked.connect(self.configure_isotopes)
        iso_layout.addWidget(self.config_iso_btn)
        
        self.iso_info = QLabel("No isotopes configured")
        self.iso_info.setStyleSheet("color: #888; padding: 5px;")
        self.iso_info.setWordWrap(True)
        iso_layout.addWidget(self.iso_info)
        
        iso_group.setLayout(iso_layout)
        left_layout.addWidget(iso_group)
        
        # 2. Detector Parameters
        det_group = QGroupBox("2. Detector Geometry")
        det_layout = QGridLayout()
        
        self.detector_params = {}
        params = [
            ('crystal_diameter', 'Crystal Ø', 60.5, 20, 100),
            ('crystal_length', 'Crystal L', 45.0, 20, 150),
            ('hole_diameter', 'Hole Ø', 9.5, 0, 30),
            ('hole_depth', 'Hole depth', 19, 0, 100),
            ('front_dead_layer', 'Dead front', 1.29, 0, 5),
            ('lateral_dead_layer', 'Dead lateral', 1.35, 0, 5),
            ('li_dead_layer', 'Li lateral', 0.8, 0, 5),
            ('li_dead_layer_front', 'Li front', 0.6, 0, 5),
            ('al_window_thickness', 'Al window', 1.4, 0, 5),
            ('al_cup_thickness', 'Al cup', 1.3, 0, 5),
            ('vacuum_gap', 'Vacuum gap', 6.0, 0, 20),
        ]
        
        for i, (key, label, default, min_v, max_v) in enumerate(params):
            row, col = i // 2, (i % 2) * 2
            det_layout.addWidget(QLabel(f"{label} (mm):"), row, col)
            spin = QDoubleSpinBox()
            spin.setRange(min_v, max_v)
            spin.setDecimals(3)
            spin.setValue(default)
            self.detector_params[key] = spin
            det_layout.addWidget(spin, row, col + 1)
        
        det_group.setLayout(det_layout)
        left_layout.addWidget(det_group)
        
        # 3. Source
        source_group = QGroupBox("3. Source (from Config tab)")
        source_layout = QVBoxLayout()
        self.source_info = QLabel("Configure in tab '1. Configuration'")
        self.source_info.setStyleSheet("color: #666;")
        source_layout.addWidget(self.source_info)
        
        refresh_btn = QPushButton("🔄 Refresh")
        refresh_btn.clicked.connect(self.load_source_config)
        source_layout.addWidget(refresh_btn)
        source_group.setLayout(source_layout)
        left_layout.addWidget(source_group)
        
        # 4. Optimization Settings
        opt_group = QGroupBox("4. Optimization Settings")
        opt_layout = QGridLayout()
        
        opt_layout.addWidget(QLabel("Parameter:"), 0, 0)
        self.param_combo = QComboBox()
        self.param_combo.addItems([
            'front_dead_layer', 'lateral_dead_layer', 
            'li_dead_layer', 'li_dead_layer_front',
            'crystal_diameter', 'crystal_length', 
            'hole_diameter', 'hole_depth',
            'al_window_thickness', 'al_cup_thickness', 'vacuum_gap'
        ])
        opt_layout.addWidget(self.param_combo, 0, 1, 1, 3)
        
        opt_layout.addWidget(QLabel("Start:"), 1, 0)
        self.start_spin = QDoubleSpinBox()
        self.start_spin.setRange(0, 100)
        self.start_spin.setDecimals(3)
        self.start_spin.setValue(0.1)
        opt_layout.addWidget(self.start_spin, 1, 1)
        
        opt_layout.addWidget(QLabel("Stop:"), 1, 2)
        self.stop_spin = QDoubleSpinBox()
        self.stop_spin.setRange(0, 100)
        self.stop_spin.setDecimals(3)
        self.stop_spin.setValue(2.0)
        opt_layout.addWidget(self.stop_spin, 1, 3)
        
        opt_layout.addWidget(QLabel("Step:"), 2, 0)
        self.step_spin = QDoubleSpinBox()
        self.step_spin.setRange(0.001, 10)
        self.step_spin.setDecimals(3)
        self.step_spin.setValue(0.2)
        opt_layout.addWidget(self.step_spin, 2, 1)
        
        opt_layout.addWidget(QLabel("Events:"), 2, 2)
        self.events_spin = QSpinBox()
        self.events_spin.setRange(1000, 10000000)
        self.events_spin.setValue(100000)
        self.events_spin.setSingleStep(10000)
        opt_layout.addWidget(self.events_spin, 2, 3)
        
        # Mode selection
        opt_layout.addWidget(QLabel("Mode:"), 3, 0)
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Mode A (cascade - use isotope)",
            "Mode B (independent - use Custom)"
        ])
        opt_layout.addWidget(self.mode_combo, 3, 1, 1, 3)
        
        self.sim_count = QLabel("Iterations: 0")
        self.sim_count.setStyleSheet("font-weight: bold; color: #E65100;")
        opt_layout.addWidget(self.sim_count, 4, 0, 1, 4)
        
        self.start_spin.valueChanged.connect(self.update_count)
        self.stop_spin.valueChanged.connect(self.update_count)
        self.step_spin.valueChanged.connect(self.update_count)
        
        opt_group.setLayout(opt_layout)
        left_layout.addWidget(opt_group)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.start_btn = QPushButton("▶️ Start Optimization")
        self.start_btn.setStyleSheet(
            "background: #4CAF50; color: white; font-weight: bold; padding: 12px;"
        )
        self.start_btn.clicked.connect(self.start_optimization)
        btn_layout.addWidget(self.start_btn)
        
        self.stop_btn = QPushButton("⏹️ Stop")
        self.stop_btn.setStyleSheet("background: #f44336; color: white; padding: 12px;")
        self.stop_btn.clicked.connect(self.stop_optimization)
        self.stop_btn.setEnabled(False)
        btn_layout.addWidget(self.stop_btn)
        
        left_layout.addLayout(btn_layout)
        
        self.progress = QProgressBar()
        self.progress.setFormat("%v/%m")
        left_layout.addWidget(self.progress)
        
        self.status = QLabel("Ready")
        left_layout.addWidget(self.status)
        
        left_layout.addStretch()
        
        # === RIGHT PANEL ===
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        tabs = QTabWidget()
        
        # Results table
        table_tab = QWidget()
        table_layout = QVBoxLayout(table_tab)
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(3)
        self.results_table.setHorizontalHeaderLabels(["Param (mm)", "χ²", "Status"])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table_layout.addWidget(self.results_table)
        
        export_btns = QHBoxLayout()
        export_csv = QPushButton("📄 Export CSV")
        export_csv.clicked.connect(self.export_csv)
        export_btns.addWidget(export_csv)
        export_best = QPushButton("⭐ Export Best")
        export_best.clicked.connect(self.export_best)
        export_btns.addWidget(export_best)
        table_layout.addLayout(export_btns)
        tabs.addTab(table_tab, "📋 Results")
        
        # Chi2 plot
        plot_tab = QWidget()
        plot_layout = QVBoxLayout(plot_tab)
        self.fig1, self.ax1 = plt.subplots(figsize=(8, 5))
        self.canvas1 = FigureCanvasQTAgg(self.fig1)
        plot_layout.addWidget(NavigationToolbar2QT(self.canvas1, self))
        plot_layout.addWidget(self.canvas1)
        tabs.addTab(plot_tab, "📈 χ² Plot")
        
        # Efficiency plot
        eff_tab = QWidget()
        eff_layout = QVBoxLayout(eff_tab)
        self.fig2, self.ax2 = plt.subplots(figsize=(8, 5))
        self.canvas2 = FigureCanvasQTAgg(self.fig2)
        eff_layout.addWidget(NavigationToolbar2QT(self.canvas2, self))
        eff_layout.addWidget(self.canvas2)
        tabs.addTab(eff_tab, "📊 Efficiency")
        
        # Log
        log_tab = QWidget()
        log_layout = QVBoxLayout(log_tab)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Courier", 9))
        log_layout.addWidget(self.log_text)
        tabs.addTab(log_tab, "📝 Log")
        
        right_layout.addWidget(tabs)
        
        # Optimal
        opt_box = QGroupBox("⭐ Optimal Configuration")
        opt_box_layout = QVBoxLayout()
        self.opt_label = QLabel("Run optimization to find optimal")
        self.opt_label.setWordWrap(True)
        opt_box_layout.addWidget(self.opt_label)
        
        self.apply_btn = QPushButton("✅ Apply Optimal")
        self.apply_btn.clicked.connect(self.apply_optimal)
        self.apply_btn.setEnabled(False)
        opt_box_layout.addWidget(self.apply_btn)
        opt_box.setLayout(opt_box_layout)
        right_layout.addWidget(opt_box)
        
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([400, 600])
        
        main_layout.addWidget(splitter)
        self.setLayout(main_layout)
        
        self.update_count()
    
    def update_count(self):
        start = self.start_spin.value()
        stop = self.stop_spin.value()
        step = self.step_spin.value()
        if step > 0:
            n = int((stop - start) / step) + 1
            n_iso = len(self.isotopes_data) if self.isotopes_data else 0
            self.sim_count.setText(f"Iterations: {n} × {n_iso} isotopes = {n * max(1, n_iso)} sims")
    
    def configure_isotopes(self):
        dialog = IsotopeEfficiencyDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            self.isotopes_data = dialog.get_data()
            
            if self.isotopes_data:
                iso_list = list(self.isotopes_data.keys())
                total_e = sum(len(v) for v in self.isotopes_data.values())
                self.iso_info.setText(
                    f"✓ {len(iso_list)} isotopes, {total_e} energies:\n" +
                    ", ".join(iso_list)
                )
                self.iso_info.setStyleSheet(
                    "color: #2E7D32; background: #E8F5E9; padding: 8px; border-radius: 4px;"
                )
            else:
                self.iso_info.setText("No isotopes configured")
                self.iso_info.setStyleSheet("color: #C62828;")
            
            self.update_count()
    
    def load_source_config(self):
        if not self.config_widget_ref:
            return
        
        cfg = self.config_widget_ref
        st = cfg.source_type.currentText()
        info = f"<b>{st}</b>"
        
        self.source_info.setText(info)
        self.source_info.setStyleSheet(
            "color: #2E7D32; background: #E8F5E9; padding: 5px; border-radius: 3px;"
        )
    
    def start_optimization(self):
        if not self.isotopes_data:
            QMessageBox.warning(self, "Error", "Configure isotopes first!")
            return
        
        if not self.sim_widget_ref:
            QMessageBox.warning(self, "Error", "Simulation widget not linked")
            return
        
        sim_path = self.sim_widget_ref.sim_path_edit.text()
        g4_path = self.sim_widget_ref.g4_path_edit.text()
        
        if not os.path.exists(sim_path):
            QMessageBox.warning(self, "Error", f"Path not found: {sim_path}")
            return
        
        # Get source config
        source_config = {}
        if self.config_widget_ref:
            cfg = self.config_widget_ref
            source_config['source_type'] = cfg.source_type.currentText()
            
            if source_config['source_type'] == 'marinelli' and cfg.marinelli_config:
                source_config.update(cfg.marinelli_config)
            elif source_config['source_type'] == 'disk' and cfg.disk_config:
                source_config.update(cfg.disk_config)
            elif source_config['source_type'] == 'volume' and cfg.cylinder_config:
                source_config.update(cfg.cylinder_config)
            else:
                source_config['source_x'] = cfg.source_x.value()
                source_config['source_y'] = cfg.source_y.value()
                source_config['source_z'] = cfg.source_z.value()
        
        mode = 'A' if 'A' in self.mode_combo.currentText() else 'B'
        
        config = {
            'param_name': self.param_combo.currentText(),
            'start': self.start_spin.value(),
            'stop': self.stop_spin.value(),
            'step': self.step_spin.value(),
            'num_events': self.events_spin.value(),
            'sim_path': sim_path,
            'geant4_path': g4_path,
            'isotopes_data': self.isotopes_data,
            'detector_params': {k: v.value() for k, v in self.detector_params.items()},
            'source_config': source_config,
            'mode': mode
        }
        
        n_params = int((config['stop'] - config['start']) / config['step']) + 1
        n_iso = len(self.isotopes_data)
        
        reply = QMessageBox.question(self, "Confirm",
            f"Parameter: {config['param_name']}\n"
            f"Range: {config['start']:.3f} → {config['stop']:.3f} mm\n"
            f"Mode: {mode}\n"
            f"Isotopes: {n_iso}\n"
            f"Total: {n_params} × {n_iso} = {n_params * n_iso} simulations\n\n"
            f"Continue?",
            QMessageBox.Yes | QMessageBox.No)
        
        if reply != QMessageBox.Yes:
            return
        
        # Reset
        self.optimization_results = []
        self.results_table.setRowCount(0)
        self.log_text.clear()
        self.progress.setMaximum(n_params)
        self.progress.setValue(0)
        
        # Start
        self.worker = OptimizationWorker(config)
        self.worker.progress.connect(self.on_progress)
        self.worker.iteration_result.connect(self.on_iteration)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.worker.log.connect(self.on_log)
        self.worker.start()
        
        self.start_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.status.setText("Running...")
        self.status.setStyleSheet("color: #FF9800; font-weight: bold;")
    
    def stop_optimization(self):
        if self.worker:
            self.worker.stop()
    
    def on_progress(self, current, total, msg):
        self.progress.setValue(current)
        self.status.setText(msg)
    
    def on_iteration(self, result):
        self.optimization_results.append(result)
        
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        
        self.results_table.setItem(row, 0, QTableWidgetItem(f"{result['param_value']:.4f}"))
        self.results_table.setItem(row, 1, QTableWidgetItem(f"{result['chi2']:.4f}"))
        
        status = QTableWidgetItem("✓")
        if result['chi2'] < 2:
            status.setBackground(QColor("#d4edda"))
        elif result['chi2'] < 10:
            status.setBackground(QColor("#fff3cd"))
        else:
            status.setBackground(QColor("#f8d7da"))
        self.results_table.setItem(row, 2, status)
        
        self.update_chi2_plot()
    
    def on_finished(self, results):
        self.start_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        
        if results.get('best_result'):
            self.best_value = results['best_value']
            self.best_param = results['param_name']
            
            self.opt_label.setText(
                f"<b>{results['param_name']}:</b> {results['best_value']:.4f} mm<br>"
                f"<b>χ²:</b> {results['best_chi2']:.4f}<br>"
                f"<b>Mode:</b> {results['mode']}"
            )
            self.opt_label.setStyleSheet("background: #E8F5E9; padding: 10px;")
            self.apply_btn.setEnabled(True)
            
            self.update_eff_plot(results['best_result'])
        
        self.status.setText("✓ Complete")
        self.status.setStyleSheet("color: #4CAF50; font-weight: bold;")
    
    def on_error(self, msg):
        self.start_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.status.setText(f"Error: {msg}")
        QMessageBox.critical(self, "Error", msg)
    
    def on_log(self, msg):
        self.log_text.append(msg)
        self.log_text.verticalScrollBar().setValue(
            self.log_text.verticalScrollBar().maximum()
        )
    
    def update_chi2_plot(self):
        if not self.optimization_results:
            return
        
        self.ax1.clear()
        
        vals = [r['param_value'] for r in self.optimization_results]
        chi2s = [r['chi2'] for r in self.optimization_results]
        
        self.ax1.plot(vals, chi2s, 'b-o', ms=8, lw=2)
        
        min_idx = np.argmin(chi2s)
        self.ax1.axvline(vals[min_idx], color='r', ls='--', lw=2)
        self.ax1.scatter([vals[min_idx]], [chi2s[min_idx]], 
                        color='red', s=200, marker='*', zorder=5)
        
        self.ax1.set_xlabel(f'{self.param_combo.currentText()} (mm)')
        self.ax1.set_ylabel('χ² reduced')
        self.ax1.set_title(f'Best: {vals[min_idx]:.4f} mm (χ²={chi2s[min_idx]:.3f})')
        self.ax1.grid(True, alpha=0.3)
        
        self.fig1.tight_layout()
        self.canvas1.draw()
    
    def update_eff_plot(self, best_result):
        self.ax2.clear()
        
        # Experimental
        for iso, peaks in self.isotopes_data.items():
            exp_e = [p['energy'] for p in peaks]
            exp_eff = [p['exp_eff'] for p in peaks]
            exp_unc = [p['exp_unc'] for p in peaks]
            self.ax2.errorbar(exp_e, exp_eff, yerr=exp_unc, fmt='o', 
                            label=f'{iso} (exp)', capsize=3)
        
        # Simulated
        if best_result['efficiencies']:
            sim_e = sorted(best_result['efficiencies'].keys())
            sim_eff = [best_result['efficiencies'][e] for e in sim_e]
            self.ax2.scatter(sim_e, sim_eff, color='black', s=100, 
                           marker='x', label='Simulated', zorder=10)
        
        self.ax2.set_xlabel('Energy (keV)')
        self.ax2.set_ylabel('Efficiency')
        self.ax2.set_title('Efficiency: Exp vs Sim')
        self.ax2.legend()
        self.ax2.grid(True, alpha=0.3)
        self.ax2.set_yscale('log')
        
        self.fig2.tight_layout()
        self.canvas2.draw()
    
    def apply_optimal(self):
        if self.best_value and self.best_param:
            self.detector_params[self.best_param].setValue(self.best_value)
            QMessageBox.information(self, "Applied", 
                                   f"{self.best_param} = {self.best_value:.4f} mm")
    
    def export_csv(self):
        if not self.optimization_results:
            return
        
        fn, _ = QFileDialog.getSaveFileName(self, "Export", "opt.csv", "CSV (*.csv)")
        if fn:
            with open(fn, 'w') as f:
                f.write("param_mm,chi2\n")
                for r in self.optimization_results:
                    f.write(f"{r['param_value']},{r['chi2']}\n")
            QMessageBox.information(self, "Exported", fn)
    
    def export_best(self):
        if not self.best_value:
            return
        
        fn, _ = QFileDialog.getSaveFileName(self, "Export", "best.txt", "Text (*.txt)")
        if fn:
            with open(fn, 'w') as f:
                f.write(f"# Optimal: {self.best_param} = {self.best_value}\n")
                for k, v in self.detector_params.items():
                    val = self.best_value if k == self.best_param else v.value()
                    f.write(f"{k} = {val}\n")
            QMessageBox.information(self, "Exported", fn)

# ==================================================================
# MAIN WINDOW
# ==================================================================

class HPGeSimulationGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HPGe Simulation Platform - Canberra Mirion Real Geometry")
        self.setGeometry(100, 100, 1400, 900)
        self.init_ui()
    
    def init_ui(self):
        cw = QWidget()
        self.setCentralWidget(cw)
        layout = QVBoxLayout(cw)
        
        self.tabs = QTabWidget()
        
        # 1. Detector Geometry (NEW - Canberra Mirion)
        self.geometry_widget = DetectorGeometryWidget()
        self.tabs.addTab(ScrollableWidget(self.geometry_widget), "1. Detector Geometry")
        
        # 2. Config
        self.config_widget = ConfigurationWidget()
        self.config_widget.set_geometry_widget(self.geometry_widget)
        self.tabs.addTab(ScrollableWidget(self.config_widget), "2. Configuration")
        
        # 3. Simulation
        self.sim_widget = SimulationControlWidget()
        self.sim_widget.set_config_widget(self.config_widget)
        self.tabs.addTab(ScrollableWidget(self.sim_widget), "3. Simulation Control")
        
        # 4. Results
        self.res_widget = EnhancedResultsWidget()
        self.res_widget.set_references(self.config_widget, self.sim_widget, self.tabs)
        self.tabs.addTab(ScrollableWidget(self.res_widget), "4. Results & COI")
        
        # 5. Geometry Optimization
        self.opt_widget = GeometryOptimizationWidget()
        self.opt_widget.set_references(self.config_widget, self.sim_widget)
        self.tabs.addTab(ScrollableWidget(self.opt_widget), "5. Geometry Optimization")
        
        layout.addWidget(self.tabs)

if __name__ == '__main__':
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    gui = HPGeSimulationGUI()
    gui.show()
    sys.exit(app.exec_())
