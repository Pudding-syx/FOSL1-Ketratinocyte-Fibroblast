## Super-Enhancer Mediated MMP3 from FOSL1+ Keratinocytes Activates Fibroblasts and Drives Immune Infiltration in Aberrant Wound Healing 

## File Structure

    \---codes_open_access_final
    +---01_scRNA_process
    |   |   environment.yaml
    |   |
    |   +---0_dnbc4tools_mapping
    |   |       dnbc4tools.txt
    |   |       reference_genome_make.txt
    |   |
    |   +---1_preprocess_integration
    |   |       1.1.LoadData.ipynb
    |   |       1.2.Preprocess.ipynb
    |   |       1.2.Preprocess.txt
    |   |       1.3.Annotation_rough.ipynb
    |   |       1.4.Figure1B-1C_S1A-S1D_S1H.ipynb
    |   |       datasheet.xlsx
    |   |
    |   +---2_keratinocyte
    |   |       2.1.keratinocyte_process.ipynb
    |   |       2.2.scenic_run.sh
    |   |       2.2.scenic_TF_analysis.ipynb
    |   |       2.3.Figure2A-I_S2A-B_S2D-G.ipynb
    |   |       keloid_KC_EMT.csv
    |   |
    |   +---3_fibroblasts
    |   |       3.1.fibroblast_process.ipynb
    |   |       3.2.Figure1D-F_S2E-G.ipynb
    |   |
    |   +---4_immune
    |   |       4.1.immune_process.ipynb
    |   |       4.2.Figure1G-J_S1I-L.ipynb
    |   |
    |   +---5_cell_cell_interaction
    |   |       5.1.cellphonedb_prepare.ipynb
    |   |       5.2.cellphone_run.ipynb
    |   |       5.2.cellphone_run_by_condition.ipynb
    |   |       5.3.cellphonedb_visualize.ipynb
    |   |       5.3.Interaction_pairs_0707.xlsx
    |   |
    |   \---myutils
    |       |   requirements.txt
    |       |   setup.py
    |       |
    |       +---build
    |       |   \---lib
    |       |       \---myutils
    |       |           |   _myutils.py
    |       |           |   __init__.py
    |       |           |
    |       |           +---plot
    |       |           |       pl.py
    |       |           |       __init__.py
    |       |           |
    |       |           +---tools
    |       |           |       tl.py
    |       |           |       __init__.py
    |       |           |
    |       |           +---_data
    |       |           |       _CellChatDB_human.py
    |       |           |       __init__.py
    |       |           |
    |       |           +---_plot
    |       |           |       pl.py
    |       |           |       _pl.py
    |       |           |       __init__.py
    |       |           |
    |       |           +---_preprocessing
    |       |           |       _doubletdetection.py
    |       |           |       _recipe.py
    |       |           |       __init__.py
    |       |           |
    |       |           \---_tools
    |       |                   _tl.py
    |       |                   __init__.py
    |       |
    |       +---myutils
    |       |   |   _myutils.py
    |       |   |   __init__.py
    |       |   |
    |       |   +---_data
    |       |   |       _CellChatDB_human.py
    |       |   |       __init__.py
    |       |   |
    |       |   +---_plot
    |       |   |       _pl.py
    |       |   |       __init__.py
    |       |   |
    |       |   +---_preprocessing
    |       |   |       _doubletdetection.py
    |       |   |       _recipe.py
    |       |   |       __init__.py
    |       |   |
    |       |   \---_tools
    |       |           _tl.py
    |       |           __init__.py
    |       |
    |       \---myutils.egg-info
    |               dependency_links.txt
    |               PKG-INFO
    |               SOURCES.txt
    |               top_level.txt
    |
    \---02_bulkRNA_process
            1.RNA-seq.ipynb
            1.string_analysis.csv
