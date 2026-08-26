#!/usr/bin/env python
# encoding: utf-8
"""

Yves NOMINE ; CBI ; 2024

06/26: added plot_all_individual_signals_pdf() to gather all per-concentration
    sensorgrammes (Data / Ref / Subtracted) in a single multi-panel PDF,
    inspired by hu_EntropyPlotAllKA.py (PdfPages + ax[r,c] grid pattern)
""" 

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# =============================================================================
# Global parameters
# =============================================================================

DOWNSAMPLING_STEP = 2
POINTS_TO_REMOVE_AT_END = 5
STEADY_STATE_THRESHOLD = 110
STEADY_STATE_WINDOW = 3

FIXED_YMIN = -1
FIXED_YMAX = 7.5
XMIN = -30.
XMAX = 500.

# --- Multi-panel PDF layout (new) -------------------------------------------
PDF_NB_ROWS = 5          # rows per page
PDF_NB_COLS = 3          # cols per page  → 15 slots/page, fits 14 sensorgrams
PDF_SUBPLOT_FONTSIZE = 7 # compact font for subplots
# ----------------------------------------------------------------------------

LIGANDS = {
    "2HPwt": {"mw": 70.0},
    "2HPnew": {"mw": 79.0},
    "2HPmut": {"mw": 70.0},
    "3AFL": {"mw": 120.0},
    "3BFL": {"mw": 119.0},
}

ANALYTES = {
    "p108a": {"mw": 0.93},
    "p108t": {"mw": 0.93},
    "p132": {"mw": 1.19},
    "p38": {"mw": 0.93},
    "p40": {"mw": 1.2},
}

DEBUG = False

plt.rcParams.update({
    'legend.fontsize': 'x-small',
    'axes.labelsize': 18,
    'axes.titlesize': 18
})


# =============================================================================
# Helper: extract concentration string from a column label
# =============================================================================
 
def _label_to_conc_str(label):
    """Return a short concentration string extracted from a column header."""
    try:
        return label.split(' ')[5].replace(',', '.')
    except IndexError:
        return label


##########################
# I/O 
##########################
def read_spr_file(filename):
    """Read and preprocess SPR data files"""
    
    file = f"{filename}_Kinetics.txt"
    data = pd.read_csv(file, sep="\t", decimal=',', encoding='ISO-8859-1') #Other possible encoding: utf-8-sig
    #data.dropna(inplace=True)
    return data.iloc[::DOWNSAMPLING_STEP].reset_index(drop=True)


def process_spr_data(sprdata, sprref, CoeffMBP, CoeffNorm):
    """
    Process SPR data and return structured results
    
    Returns:
    list of dicts containing processed data for each column pair
    """
    
    num_pairs = int(sprdata.shape[1]/2)
    colData_headers = list(sprdata.columns)
    processed_data = []
    print(f"Number of pairs: {num_pairs}")
    
    for i in range(num_pairs):
        end_index = -POINTS_TO_REMOVE_AT_END
        data = {
            'x_values': sprdata.iloc[:end_index, 2*i],
            'y_Data': sprdata.iloc[:end_index, 2*i + 1],
            #'y_Ref': sprref.iloc[:end_index, 2*i + 1], #See in loop if len(sprref)
            'labelData': colData_headers[2*i + 1],
            #'labelRef' : colRef_headers[2*i + 1]
        }
        # print(f"Index in Num_Pairs loop: {i} with CoeffMBP = {CoeffMBP} and CoeffNorm = {CoeffNorm}")

        if len(sprref) is not None:
            colRef_headers  = list(sprref.columns)
            data['labelRef'] = colRef_headers[2*i + 1]
            data['y_Ref'] = sprref.iloc[:end_index, 2*i + 1]
            data['y_subtracted'] = ( data['y_Data'] - data['y_Ref']*CoeffMBP )
        else:
            data['y_subtracted'] = data['y_Data'] 
        data['y_subtracted_Norm'] = data['y_subtracted'] / CoeffNorm
        processed_data.append(data)
        #print(f"Processed data: {processed_data}")
    
        concentrationData = float(data['labelData'].split(' ')[5].replace(',', '.'))
        if len(sprref) is not None:
            concentrationRef = float(data['labelRef'].split(' ')[5].replace(',', '.'))
            if concentrationData != concentrationRef:
                print(f"Concentration for Data {concentrationData} vs. Concentration for Ref {concentrationRef}")
                raise ValueError("Data to be processed are coming with different concentrations. Check the input files.")
        
    return processed_data


def calculate_averaged_signal(data, threshold=STEADY_STATE_THRESHOLD):
    """
    Calculate the averaged signal around the first x-value above a threshold.
    
    Parameters:
        data (dict): Processed dataset with 'x_values' and 'y_subtracted_Norm'.
        threshold (float): Threshold value for x_values.
    
    Returns:
        float, float: Concentration and averaged signal.
    """
    x_values = np.asarray(data['x_values'])
    y_values = data['y_subtracted_Norm']
    
    indices = np.where(x_values > threshold)[0]
    if len(indices) == 0:
        raise ValueError(
            f"No point found above threshold {threshold}"
        )

    ind = indices[0]
    
    if DEBUG:
        print(ind)
        print(x_values[ind-10:ind+10])
        print(y_values.iloc[ind-2:ind+3])
    
    start = max(0, ind - STEADY_STATE_WINDOW)
    stop  = ind + STEADY_STATE_WINDOW + 1
    avg_signal = y_values.iloc[start:stop].mean()
    
    concentration = float(data['labelData'].split(' ')[5].replace(',', '.'))
    
    return concentration, avg_signal


def plot_individual_signal(data, ref, index):
    """
    Plot individual SPR signal with subtracted reference.
    
    Parameters:
        data (dict): Processed dataset
        ref: 1 to plot the ref, 0 if not
        index (int): Index of the dataset for plot labeling.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(data['x_values'], data['y_Data'], 'g', label="Data")
    if ref == 1:
        plt.plot(data['x_values'], data['y_Ref'], 'b', label="RefMBP")
    plt.plot(data['x_values'], data['y_subtracted'], 'r', label="Subtracted")
    plt.xlabel('Time (s)')
    plt.ylabel('Response Units')
    plt.title(f"Signal {index}: {data['labelData']}")
    plt.legend(loc='upper right')
    plt.grid(True)
    #FigName = "SteadyState_" + data['labelData'] + ".eps"
    #plt.savefig(FigName)
    plt.show()


def plot_final_signals(processed_spr, title, FigName):
    """
    Plot all subtracted / normalized signals superimposed on a single axes.
    
    Parameters:
        processed_spr (list[dict]): List of processed datasets.
        title (str): Plot title.
        FigName (str): Output file path
    """
    
    plt.figure(figsize=(12, 8))
    cmap = plt.get_cmap('jet')  
    colors = [cmap(i) for i in np.linspace(0, 1, len(processed_spr))]
    max_value = 0
    min_value = 0
    
    for idx, data in enumerate(processed_spr):
        if DEBUG:
            print(f"LenX: {len(data['x_values'])},  LenDataY: {len(data['y_Data'])},  LenY_DataSub: {len(data['y_subtracted'])},  LenY_Ref: {len(data['y_Ref'])} and LenY_Norm: {len(data['y_subtracted_Norm'])}")
        plt.plot(data['x_values'], data['y_subtracted_Norm'], label=data['labelData'], color=colors[idx])
        max_value = max(max_value, data['y_subtracted_Norm'].max())    
        min_value = min(min_value, data['y_subtracted_Norm'].min())    

    Upper = FIXED_YMAX
    Lower = FIXED_YMIN
    if DEBUG:
        print(Lower, Upper)

    plt.xlabel('Time (s)')
    plt.ylabel('Normalized Response Units')
    plt.legend(loc='upper right')
    plt.axis([XMIN, XMAX, Lower, Upper])
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(FigName)
    plt.show()


# =============================================================================
# Plotting — multi-panel PDF (one PDF for all concentrations)
# =============================================================================
 
def plot_all_individual_signals_pdf(processed_spr, ref, title, pdf_filename,NbRow=PDF_NB_ROWS, NbCol=PDF_NB_COLS):
    """
    Save all per-concentration sensorgrammes (Data / Ref / Subtracted) into a
    single multi-panel PDF, using a NbRow × NbCol subplot grid per page.
 
    Parameters:
        processed_spr (list[dict]): List of processed datasets from process_spr_data().
        ref (int):        1 to overlay the reference channel, 0 otherwise.
        title (str):      Overall figure suptitle (dataset identifier).
        pdf_filename (str): Output PDF path.
        NbRow (int):      Number of subplot rows per page  (default 5).
        NbCol (int):      Number of subplot columns per page (default 3).
    """
    fls       = PDF_SUBPLOT_FONTSIZE
    slots     = NbRow * NbCol          # panels available per page
    page_idx  = 1                      # for potential debug / filenames
    r, c      = 0, 0                   # current row / col position
    fig, axes = None, None             # created lazily: only when a plot is needed
 
    with PdfPages(pdf_filename) as pdf:
        for idx, data in enumerate(processed_spr):
            
            # ── open a new page only when actually needed ────────────────────
            if fig is None:
                fig, axes = plt.subplots(NbRow, NbCol,
                                         figsize=(NbCol * 5, NbRow * 3),
                                         squeeze=False)   # always 2-D array
                fig.suptitle(title, fontsize=fls + 2, fontweight='bold')
                
            ax = axes[r, c]
 
            # ── traces ──────────────────────────────────────────────────────
            ax.plot(data['x_values'], data['y_Data'],
                    'g', linewidth=0.8, label="Data")
            if ref == 1 and 'y_Ref' in data:
                ax.plot(data['x_values'], data['y_Ref'], 'b', linewidth=0.8, label="RefMBP")
            ax.plot(data['x_values'], data['y_subtracted'], 'r', linewidth=0.8, label="Subtracted")
 
            # ── cosmetics ───────────────────────────────────────────────────
            conc_str = _label_to_conc_str(data['labelData'])
            ax.set_title(f"#{idx + 1}  [{conc_str} uM]", fontsize=fls, fontweight='bold')
            ax.set_xlabel('Time (s)', fontsize=fls)
            if c == 0:
                ax.set_ylabel('Response Units', fontsize=fls)
            ax.tick_params(labelsize=fls - 1)
            ax.legend(fontsize=fls - 1, loc='upper right')
            ax.grid(True, linewidth=0.4)
 
            if DEBUG:
                print(f"  page {page_idx}, subplot [{r},{c}] → {data['labelData']}")
 
            # ── advance position ─────────────────────────────────────────────
            c += 1
            if c == NbCol:
                c = 0
                r += 1
                if r == NbRow:
                    # Current page is full → save it and open a new one
                    fig.tight_layout(rect=[0, 0, 1, 0.96])
                    pdf.savefig(fig)
                    plt.close(fig)
                    page_idx += 1
                    r, c = 0, 0
                    fig, axes = None, None
 
        # ── last (possibly partial) page ─────────────────────────────────────
        # If r==0 and c==0 the last page was already saved full (exactly
        # NbRow*NbCol plots): just discard the empty figure that was opened
        # speculatively after that save.
        if fig is not None:
            used_on_last_page = idx % slots + 1
            for empty in range(used_on_last_page, slots):
                axes[empty // NbCol, empty % NbCol].set_visible(False)
        
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            pdf.savefig(fig)
            plt.close(fig)
 
    print(f"Multi-panel PDF saved → {pdf_filename}")
 

def MAIN_PROCESS(file1, file2, CoeffMBP=1.0, CoeffNorm = 1.0, Immob = 1000., plot_individual_pdf=True):
    """
    Main function to process and plot SPR data.
    
    Parameters:
        file1 (str): First data file name (without extension).
        file2 (str): Second data file name (without extension).
        CoeffMBP (float): MBP reference correction coefficient.
        CoeffNorm (float): Normalization coefficient.
        Immob (float): Ligand immobilisation level in RU; used to compute CoeffNorm_A_L = Immob × Analyte_MW / Ligand_MW.
        plot_individual_pdf (bool): When True (default), generate the multi-panel PDF with all per-concentration sensorgrams.
            Set to False to skip it.
    """
    
    if DEBUG:
        print(f"File1: {file1} and File2: {file2}")
        print(f"Len of File1: {len(file1)} and Length of File2: {len(file2)}")
        
    # Read data files
    SPRdata = read_spr_file(file1)
    if len(file2) != 0:
        SPRref = read_spr_file(file2)
        ref = 1
    else:
        SPRref = None
        ref = 0


    # ––––––––– Calculate the Normalization Coefficient (based on analyte MW, ligand MW and ligand immob level) –––––
    Chain1        = file1.split('_')
    Chain2        = file2.split('_')
    Analyte       = Chain1[2]     #ex : p132
    Ligand        = Chain1[4]      #ex : 3BFL
    Analyte_MW    = ANALYTES[Analyte]["mw"]
    Ligand_MW     = LIGANDS[Ligand]["mw"]
    CoeffNorm_A_L = Immob * Analyte_MW / Ligand_MW
    CoeffNorm     = CoeffNorm_A_L
    print(f"Analyte: {Analyte} ({Analyte_MW} kDa ); Ligand: {Ligand} ({Ligand_MW} kDa)")
    
    
    #if SPRdata.shape != SPRref.shape:
    #    raise ValueError("Files have different structures. Check the input files.")


    # ––––––––– Process SPR data ––––––––– 
    processed_spr = process_spr_data(SPRdata, SPRref, CoeffMBP, CoeffNorm)
    if DEBUG:
        print(processed_spr)


    # ––––––––– Steady-state file ––––––––– 
    with open(f"{file1}_SteadyState_NormAL.txt", "w") as target:
        target.write("Concentration  AveragedSignal\n")
        
        for index, data in enumerate(processed_spr):
            try:
                # Calculate averaged signal
                Conc, avg_signal = calculate_averaged_signal(data)
                
                # Save to file
                target.write(f"{Conc}\t{avg_signal}\n")
                print(f"Concentration: {Conc}, Averaged Signal: {avg_signal}")
                
                # Plot individual signal
                # plot_individual_signal(data, ref, index + 1)
            
            except ValueError as e:
                print(f"Skipping dataset {index + 1}: {e}")
    
    #  ––––––––– Build title string  ––––––––– 
    if len(file2) != 0:
        TITLE = Chain1[0] + '_' + Chain1[1] + '_' + Chain1[2] + '_' + Chain1[4] + '-' + Chain2[4] + '_' + "CoeffMBP" + str(CoeffMBP) + '_' + "CoeffNorm_A_L" + str(CoeffNorm_A_L)
    else:
        TITLE = Chain1[0] + '_' + Chain1[1] + '_' + Chain1[2] + '_' + Chain1[4] + '_' + "CoeffMBP" + str(CoeffMBP) + '_' + "CoeffNorm_A_L" + str(CoeffNorm_A_L)
    
    print(Chain1, Chain2, TITLE)


    # ── multi-panel PDF: one page with all per-concentration sensorgrammes ────
    if plot_individual_pdf:
        figPDF_individual = (f"{file1}_Kinetics_Individual_MBP.pdf")
        plot_all_individual_signals_pdf(processed_spr, ref, TITLE, figPDF_individual, NbRow=PDF_NB_ROWS, NbCol=PDF_NB_COLS)


    # ––––––––– Superimposed normalized signals ––––––––– 
    #figPNG = file1 + '_Kinetics_MBPsubtracted_Improved_SameScale.png'
    figEPS = file1 + '_Kinetics_MBPsubtracted_Improved_SameScale_NormAL.eps'
    
    plot_final_signals(processed_spr, TITLE, figEPS)

# ————————— Entry Point —————————————————

if __name__ == "__main__":
    ###
    # CoeffMBP and CoeffNorm were determined separatly
    ###
    
    MAIN_PROCESS("241213_D_p108a_Fc4-1corr_3AFL", "241213_D_p108a_Fc2-1corr_MBP", CoeffMBP = 0.35, Immob = 1050)


