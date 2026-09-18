# Transcriptomics Notebook

**Course:** Into to Ecological Genomics Fall 2026

**Name:** Kai Roberts-Speaker

------------------------------------------------------------------------

## 9/15/2026 - Setting up Lab Notebook and Setting up Markdown

-   Setting up transcriptomics notebook

-   Learn to make notes in markdown

-   Push notes to github

    **Working Directory**

    `/gpfs1/home/k/i/kirobert/Projects/eco_genomics_2026/transcriptomics`

    **Input Files**

    `none`

    **Output Files**

    `gpfs1/home/k/i/kirobert/Projects/eco_genomics_2026/transcriptomics/TRANSCRIPTOMICS_NOTEBOOK.md`

    **Programs and Dependenceis**

-   `R Version 4.5.1`

-   `R Studio`

    **Scripts**

-   `none`

    **Code**

    `` the usage of ##_before words determines what size font will be produced, where more hashtags imply  a smaller size. **_** around words will bold them. *** will create a line. /table will create a table /image will allow import of an image. `_` will create a textbowx ``

**Table**

| Col1 | Col2 | Col3                    |
|------|------|-------------------------|
|      |      |                         |
|      |      |                         |
|      |      | ![](images/images.webp) |

------------------------------------------------------------------------

# Transcriptomics Notebook

**Course:** Into to Ecological Genomics Fall 2026

**Name:** Kai Roberts-Speaker

------------------------------------------------------------------------

## 9/17/2026 Learning to use Commands on the Vacc Linux Server

-   Set path to enter the "cl" directory, before navigating to the "biol3990" and subsequently "Transcriptomics" folders and opening the "CleanData" directory which contains cleaned RNAseq data files in fq.gz form

-   Use commands ll and ls to list the detailed and simplified contents of a directory

-   use command "zcat ..." (where ... is a .gz) file to unzip and \| to send output to command "head -n 4" to read the first 4 lines of file contents

-   **Working Directory**

    `/gpfs1/cl/biol3990/Transcriptomics/CleanData`

    **Input Files**

    `AA_F0_Rep1_2_clean.fq.gz`

    **Output Files**

    `none`

    **Programs and Dependenceis**

-   `VACC Linux Server`

-   `Bash`

    **Scripts**

-   `none`

    **Code**

    ```         
    cd /gpfs1/cl/biol3990 to set working directory to biol3990

    cd Transcriptomics to enter transcriptomics

    ll to view the contents of a directory and info about said directory incluiding ownership and size

    ls to list contents of the directory 

    command zcat ...gz to uncompress a .gz file

    ... | ... to send the command to the following function (pipe)

     "head" begins reading from the start of the file

    -n 4 to quantify reading the first 4 lines of the file

    zcat...gz | wc -l  where the data in the  unzipped file is analyzed  by the "word count" command which is modified as wc -l where -l implies quantifying the number of lines. Output represents 4x number of reads as each read contains 4 lines
    ```

------------------------------------------------------------------------
