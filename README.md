# 🌐 ChemInsight

## Overview
ChemInsight is a module of **[AIDrugApp v1.2.6](https://aidrugapp.streamlit.app/)** that facilitates molecule identification by:
- Converting chemical names to SMILES
- Identifying compound names and their 2D structures from SMILES strings
- Detecting molecular similarities from user data

## 🔧 Features
- **Single or Batch Processing**: Search for a single molecule or process multiple molecules via batch file upload.
- **Flexible Identification**:
  - Convert **Name to SMILES** (Canonical or Isomeric)
  - Convert **SMILES to Compound Information**
  - **Molecular Similarity** detection
  - **2D Structure Rendering**
- **CSV Support**: Upload a CSV file for batch processing.
- **Downloadable Outputs**: Save results in CSV format.

## 📌 How to Use ChemInsight

### Step 1: Choose Search Type
Select whether to process a **single molecule** or **multiple molecules (batch processing)**.

### Step 2: Select Identification Type
Choose the type of data retrieval:
- **Name to SMILES** (Canonical or Isomeric)
- **SMILES to Compound Info**
- **Molecular Similarity**
- **2D Structure Rendering**

### Step 3: Provide Input Data
- **Single Molecule**: Enter the molecular name or SMILES string.
- **Batch Processing**: Upload a CSV file (Example input file available [here](https://github.com/DivyaKarade/Example-.csv-input-files--AIDrugApp-v1.2/blob/main/Mol_identifier.csv)).

### Step 4: Retrieve Results
Click the corresponding button to generate results, which can be copied or downloaded as a CSV file.

## 🛠 Dependencies
This application requires the following Python libraries:
```bash
pip install matplotlib pandas streamlit rdkit pubchempy mols2grid
```

## 🚀 Running the Application
To run ChemInsight locally, use:
```bash
streamlit run chem_insight.py
```

## 📄 Reference
- **Publication**: [Custom ML Module of AIDrugApp](https://doi.org/10.33774/chemrxiv-2021-3f1f9) by Divya Karade

## 🖥️ Screenshot
![ChemInsight Screenshot](https://your-image-link-here.com/screenshot.png)
