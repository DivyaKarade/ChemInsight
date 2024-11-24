"""
  Copyright 2024 Divya Karade

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from io import BytesIO
import streamlit.components.v1 as components
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, AllChem, MACCSkeys
import pubchempy as pcp
import mols2grid
import urllib

# Page expands to full width
st.set_page_config(page_title='AIDrugApp', page_icon='🌐', layout="wide")

st.title('ChemInsight')
st.success(
    "ChemInsight is a module of [**AIDrugApp v1.2.6**](https://aidrugapp.streamlit.app/) that facilitates molecule "
    "identification by converting chemical names to SMILES, identifying compound names and their 2D structures from "
    "SMILES strings, and detecting molecular similarities from user data.")

expander_bar = st.expander("👉 More information")
expander_bar.markdown("""
            * **Python libraries:** pubchempy, pandas, rdkit, mols2grid, matplotlib
            * **Publications:** Divya Karade. (2021). Custom ML Module of AIDrugApp for Molecular Identification, Descriptor Calculation, and Building ML/DL QSAR Models. [ChemRxiv Preprint] (https://doi.org/10.33774/chemrxiv-2021-3f1f9).
            """)

expander_bar = st.expander("👉 How to use ChemInsight?")
expander_bar.markdown("""
                **Step 1:** On the "User Input Panel" first select whether you would like to search for a single molecule or upload a batch file for multiple molecules (Example input batch file given)
                """)
expander_bar.markdown("""
                **Step 2:** Select type of data to identify
                """)
expander_bar.markdown("""
                **Step 3:** Input molecular names or SMILES as per requirement in single string or batch file
                    """)
expander_bar.markdown("""
                **Step 4:** Click on the "button" and the results will be displayed below to copy or download
                """)

# Sidebar
# Collects user input features into dataframe
st.sidebar.header('⚙️ USER INPUT PANEL')
st.sidebar.subheader('1. Type of Input data')
add_selectbox1 = st.sidebar.radio(
    "How would you like to search?",
    ("Single molecule", "Multiple molecules (Batch)"))

st.sidebar.subheader('2. Type of retrieving data')
add_selectbox2 = st.sidebar.radio(
    "What would you like to identify?",
    ("Name to SMILE", "SMILE to Compound", "Molecular Similarity", "2-D Structure"))

if add_selectbox2 == "Name to SMILE":
    add_selectbox3 = st.sidebar.selectbox("Select type of SMILES to retrieve",
                                          ("Canonical SMILES", "Isomeric SMILES"))

if add_selectbox1 == 'Single molecule':
    if add_selectbox2 == "Name to SMILE":
        # Get the feature input from the user
        st.sidebar.subheader('3. Enter compound name')
        # def user_input_features():
        mol_name = st.sidebar.text_input('Enter compound name to retrieve SMILES', 'compound name')

        # Store a dictionary into a variable
        user_data = {'Name': mol_name}

        # Transform a data into a dataframe
        user_input = pd.DataFrame(user_data, index=[0])
        df = pd.concat([user_input['Name']], axis=1)
        df.to_csv('molecule.txt', sep='\t', header=False, index=False)
        st.write(df)
        List_of_Chemicals = df['Name'].values

        if st.sidebar.button("😊 GET SMILES"):
            if add_selectbox3 == "Canonical SMILES":
                # list of chemical names
                for chemical_name in List_of_Chemicals:
                    cid = pcp.get_cids(chemical_name)
                    prop2 = pcp.get_compounds(chemical_name, 'name')
                    for compound in prop2:
                        x = compound.canonical_smiles
                        x1 = (chemical_name + ' ' + str(x))
                        st.subheader('Canonical SMILES')
                        st.write(x1)

            if add_selectbox3 == "Isomeric SMILES":
                # list of chemical names
                for chemical_name in List_of_Chemicals:
                    cid = pcp.get_cids(chemical_name)
                    prop2 = pcp.get_compounds(chemical_name, 'name')
                    for compound in prop2:
                        x = compound.isomeric_smiles
                        x1 = (chemical_name + ' ' + str(x))
                        st.subheader('Isomeric SMILES')
                        st.write(x1)

    if add_selectbox2 == "SMILE to Compound":
        # Get the feature input from the user
        st.sidebar.subheader('3. Enter Canonical SMILES')
        # def user_input_features():
        smiles = st.sidebar.text_input('Enter SMILES to retrieve compound info', 'Canonical SMILES')

        # Store a dictionary into a variable
        user_data = {'smiles': smiles}

        # Transform a data into a dataframe
        user_input = pd.DataFrame(user_data, index=[0])
        df = pd.concat([user_input['smiles']], axis=1)
        df.to_csv('molecule.smi', sep='\t', header=False, index=False)
        st.write(df)
        smile = df['smiles'].values

        if st.sidebar.button("😊 GET COMPOUND"):
            # list of chemical names
            for smiles in smile:
                prop2 = pcp.get_compounds(smiles, 'smiles')
                comp = (smiles + ' ' + str(prop2))
                st.subheader('PubChem Compound ID')
                st.write(comp)
                for compound in prop2:
                    df4 = pcp.compounds_to_frame(prop2,
                                                 properties=['synonyms', 'canonical_smiles',
                                                             'molecular_formula',
                                                             'iupac_name', 'inchi', 'inchikey'])
                    st.subheader('Compound Info')
                    st.write(df4)
                    st.download_button('Download CSV', df4.to_csv(), 'Cmpd_info.csv', 'text/csv')

    if add_selectbox2 == "Molecular Similarity":
        # Get the feature input from the user
        st.sidebar.subheader('3. Enter SMILES to get molecular similarity')
        # def user_input_features():
        smiles1 = st.sidebar.text_input('Enter SMILES1', 'Canonical SMILES1')
        smiles2 = st.sidebar.text_input('Enter SMILES2', 'Canonical SMILES2')

        # Store a dictionary into a variable
        user_data1 = {'smiles1': smiles1}
        user_data2 = {'smiles2': smiles2}

        # Transform a data into a dataframe
        user_input1 = pd.DataFrame(user_data1, index=[0])
        user_input2 = pd.DataFrame(user_data2, index=[0])
        df = pd.concat([user_input1['smiles1'], user_input2['smiles2']], axis=1)
        # df2 = pd.concat([df['smiles']], axis=1)
        df.to_csv('molecule.smi', sep='\t', header=False, index=False)
        st.write(df)
        # s1 = df['smiles1'].values
        # s2 = df['smiles2'].values

        if st.sidebar.button("😊 GET SIMILARITY"):
            # list
            st.subheader('Similarity Scores')

            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)

            st.write("Tanimoto    :", round(DataStructs.TanimotoSimilarity(fp1, fp2), 4))
            st.write("Dice        :", round(DataStructs.DiceSimilarity(fp1, fp2), 4))
            st.write("Cosine      :", round(DataStructs.CosineSimilarity(fp1, fp2), 4))
            st.write("Sokal       :", round(DataStructs.SokalSimilarity(fp1, fp2), 4))
            st.write("McConnaughey:", round(DataStructs.McConnaugheySimilarity(fp1, fp2), 4))

    if add_selectbox2 == "2-D Structure":
        # Get the feature input from the user
        st.sidebar.subheader('3. Enter Canonical SMILES')
        # def user_input_features():
        smiles = st.sidebar.text_input('Enter SMILES to get 2-D structure', 'Canonical SMILES')

        # Store a dictionary into a variable
        user_data = {'smiles': smiles}

        # Transform a data into a dataframe
        user_input = pd.DataFrame(user_data, index=[0])
        df = pd.concat([user_input['smiles']], axis=1)
        df.to_csv('molecule.smi', sep='\t', header=False, index=False)
        st.write(df)
        smile = df['smiles'].values

        if st.sidebar.button("😊 GET STRUCTURE"):
            st.subheader('2-D structure')
            mol = Chem.MolFromSmiles(smiles)
            # x = st.image(Draw.MolToImage(mol))
            img = (Draw.MolToImage(mol))
            bio = BytesIO()
            img.save(bio, format='png')
            st.image(img)

if add_selectbox1 == 'Multiple molecules (Batch)':
    if add_selectbox2 == "Name to SMILE":
        # Sidebar
        with st.sidebar.subheader('3. Upload your CSV data'):
            uploaded_file = st.sidebar.file_uploader(
                "Upload your input CSV file containing 'Name' as one column of compound names", type=["csv"])
            st.sidebar.markdown("""
        [Example CSV input file](https://github.com/DivyaKarade/Example-.csv-input-files--AIDrugApp-v1.2/blob/main/Mol_identifier.csv)
        """)

        if uploaded_file is not None:
            # Read CSV data
            @st.cache_data
            def load_csv():
                csv = pd.read_csv(uploaded_file)
                return csv


            df1 = load_csv()
            df2 = df1.iloc[:, 0:8]
            X = df2
            # Write CSV data
            df2.to_csv('molecule.txt', sep='\t', header=False, index=False)
            st.subheader('Uploaded data')
            st.write(df2)
            List_of_Chemicals = df2['Name'].values

            if st.sidebar.button("😊 GET SMILES"):
                if add_selectbox3 == "Canonical SMILES":
                    st.subheader('Canonical SMILES')
                    # list of chemical names
                    data = []
                    for Name in List_of_Chemicals:
                        try:
                            df4 = pcp.get_properties(
                                ['canonical_smiles'], Name, 'name')

                            # st.write(df4)
                            data.append(df4)
                        except (pcp.BadRequestError, TimeoutError, urllib.error.URLError, ValueError):
                            pass
                    # st.write(data)

                    rows = []
                    columns = data[0][0].keys()
                    for i in range(len(data)):
                        rows.append(data[i][0].values())
                    props_df = pd.DataFrame(data=rows, columns=columns)
                    # st.write(props_df)
                    props_df.insert(0, 'Name', df2['Name'], True)
                    st.write(props_df)
                    st.download_button('Download CSV', props_df.to_csv(), 'Can_smile.csv', 'text/csv')

                if add_selectbox3 == "Isomeric SMILES":
                    # list of chemical names
                    st.subheader('Isomeric SMILES')
                    # list of chemical names
                    data = []
                    for Name in List_of_Chemicals:
                        try:
                            df4 = pcp.get_properties(
                                ['isomeric_smiles'], Name, 'name')

                            # st.write(df4)
                            data.append(df4)
                        except (pcp.BadRequestError, TimeoutError, urllib.error.URLError, ValueError):
                            pass
                    # st.write(data)

                    rows = []
                    columns = data[0][0].keys()
                    for i in range(len(data)):
                        rows.append(data[i][0].values())
                    props_df = pd.DataFrame(data=rows, columns=columns)
                    # st.write(props_df)
                    props_df.insert(0, 'Name', df2['Name'], True)
                    st.write(props_df)
                    st.download_button('Download CSV', props_df.to_csv(), 'isomeric_smile.csv', 'text/csv')

    if add_selectbox2 == "SMILE to Compound":
        # Sidebar
        with st.sidebar.subheader('3. Upload your CSV data'):
            uploaded_file = st.sidebar.file_uploader(
                "Upload your input CSV file containing 'smile' as one column of molecular SMILES", type=["csv"])
            st.sidebar.markdown("""
            [Example CSV input file](https://github.com/DivyaKarade/Example-.csv-input-files--AIDrugApp-v1.2/blob/main/Mol_identifier.csv)
            """)

        if uploaded_file is not None:
            # Read CSV data
            @st.cache_data
            def load_csv():
                csv = pd.read_csv(uploaded_file)
                return csv


            df1 = load_csv()
            df2 = df1.iloc[:, 0:]
            X = df2
            # Write CSV data
            df2.to_csv('molecule.smi', sep='\t', header=False, index=False)
            st.subheader('Uploaded data')
            st.write(df2)
            smile = df2['smiles'].values

            if st.sidebar.button("😊 GET COMPOUND"):

                # list of chemical names
                st.subheader('PubChem Compound Information')

                data = []
                for smiles in smile:
                    try:
                        df4 = pcp.get_properties(
                            ['canonical_smiles', 'molecular_formula', 'iupac_name', 'inchi', 'inchikey'],
                            smiles,
                            'smiles')
                        # st.write(df4)
                        data.append(df4)
                    except (pcp.BadRequestError, TimeoutError, urllib.error.URLError, ValueError):
                        pass
                # st.write(data)

                rows = []
                columns = data[0][0].keys()
                for i in range(len(data)):
                    rows.append(data[i][0].values())
                props_df = pd.DataFrame(data=rows, columns=columns)
                # st.write(props_df)
                props_df.insert(0, 'smiles', df2['smiles'], True)
                st.write(props_df)
                st.download_button('Download CSV', props_df.to_csv(), 'Cmpd_info.csv', 'text/csv')

    if add_selectbox2 == "Molecular Similarity":
        # Sidebar
        with st.sidebar.subheader('3. Upload your CSV data'):
            uploaded_file = st.sidebar.file_uploader(
                "Upload input file containing 'Name' & 'smile' as two columns of compound names and SMILES",
                type=["csv"])
            st.sidebar.markdown("""
                            [Example CSV input file](https://github.com/DivyaKarade/Example-.csv-input-files--AIDrugApp-v1.2/blob/main/Mol_identifier.csv)
                            """)

        if uploaded_file is not None:
            # Read CSV data
            @st.cache_data
            def load_csv():
                csv = pd.read_csv(uploaded_file)
                return csv


            df1 = load_csv()
            df2 = df1.iloc[:, 0:]
            X = df2
            # Write CSV data
            df2.to_csv('molecule.smi', sep='\t', header=False, index=False)
            st.subheader('Uploaded data')
            st.write(df2)
            smile = df2['smiles'].values
            Name = df2['Name'].values

            if st.sidebar.button("😊 GET SIMILARITY"):
                st.write("****")
                st.subheader('Computation of similarity scores')
                st.info('**Similarity scores between compounds**')

                mols = [Chem.MolFromSmiles(smiles) for smiles in smile]
                fps = [MACCSkeys.GenMACCSKeys(x) for x in mols]
                st.write("Number of compounds:", len(mols))
                st.write("Number of fingerprints:", len(fps))
                st.write("The number of compound pairs:", (len(fps) * (len(fps) - 1)) / 2)

                scores = []

                for i in range(0, len(fps)):

                    if i == 0:
                        print("Processing compound ", end='')

                    if i % 100 == 0:
                        print(i, end=' ')

                    for j in range(i + 1, len(fps)):
                        scores.append(DataStructs.FingerprintSimilarity(fps[i], fps[j]))

                st.write("Number of scores : ", len(scores))
                st.write("****")
                st.info('**Similarity Scores**')
                st.write("Tanimoto    :", round(DataStructs.TanimotoSimilarity(fps[0], fps[1]), 4))
                st.write("Dice        :", round(DataStructs.DiceSimilarity(fps[0], fps[1]), 4))
                st.write("Cosine      :", round(DataStructs.CosineSimilarity(fps[0], fps[1]), 4))
                st.write("Sokal       :", round(DataStructs.SokalSimilarity(fps[0], fps[1]), 4))
                st.write("McConnaughey:", round(DataStructs.McConnaugheySimilarity(fps[0], fps[1]), 4))

                # Generate a histogram that shows the distribution of the pair-wise scores
                st.write("****")
                st.subheader('Distribution of similarity scores')
                st.info('**Histograms showing the distribution of the pair-wise scores**')
                mybins = [x * 0.01 for x in range(101)]

                # Create the figure and axes
                fig, ax = plt.subplots(figsize=(8, 4), dpi=300)

                # First subplot: Distribution
                ax.subplot(1, 2, 1)
                ax.set_title("Distribution")
                ax.hist(scores, bins=mybins)

                # Second subplot: Cumulative Distribution
                ax.subplot(1, 2, 2)
                ax.set_title("Cumulative Distribution")
                ax.hist(scores, bins=mybins, density=True, cumulative=1)
                ax.plot([0, 1], [0.95, 0.95])

                # Pass the figure to st.pyplot()
                st.pyplot(fig)

                st.info("**Interpretation of similarity scores**")
                st.write(
                    "Average similarity score between two compounds (computed using the Tanimoto equation and MACCS keys) :",
                    sum(scores) / len(scores))
                # to find a threshold for top 3% compound pairs (i.e., 97% percentile)
                st.write("Total compound pairs:   ", len(scores))
                st.write("95% of compound pairs:  ", len(scores) * 0.97)
                st.write("Score at 95% percentile:", scores[round(len(scores) * 0.97)])

                st.info(
                    '**Table showing the % of compound pairs and their similarity scores.** (Table Columns: Similarity score, Number of compound pairs, % of compound pairs)')
                for i in range(21):
                    thresh = i / 20
                    num_similar_pairs = len([x for x in scores if x >= thresh])
                    prob = num_similar_pairs / len(scores) * 100
                    st.write("%.3f  %8d  (%8.4f %%)" % (thresh, num_similar_pairs, round(prob, 4)))

                st.write("****")
                st.info('**Pair-wise similarity scores among molecules**'
                        ' (To make higher scores easier to find, they are indicated with the "+" character(s).)')
                for i in range(0, len(fps)):
                    for j in range(i + 1, len(fps)):

                        score = DataStructs.FingerprintSimilarity(fps[i], fps[j])
                        st.write(Name[i], "vs.", Name[j], ":", round(score, 3), end='')

                        if score >= 0.85:
                            st.write(" ++++ ")
                        elif score >= 0.75:
                            st.write(" +++ ")
                        elif score >= 0.65:
                            st.write(" ++ ")
                        elif score >= 0.55:
                            st.write(" + ")
                        else:
                            st.write(" ")

    if add_selectbox2 == "2-D Structure":
        # Sidebar
        with st.sidebar.subheader('3. Upload your CSV data'):
            uploaded_file = st.sidebar.file_uploader(
                "Upload your input CSV file containing 'smile' as one column of molecular SMILES", type=["csv"])
            st.sidebar.markdown("""
                    [Example CSV input file](https://github.com/DivyaKarade/Example-.csv-input-files--AIDrugApp-v1.2/blob/main/Mol_identifier.csv)
                    """)

        if uploaded_file is not None:
            # Read CSV data
            @st.cache_data
            def load_csv():
                csv = pd.read_csv(uploaded_file)
                return csv


            df1 = load_csv()
            df2 = df1.iloc[:, 0:]
            X = df2
            # Write CSV data
            df2.to_csv('molecule.smi', sep='\t', header=False, index=False)
            st.subheader('Uploaded data')
            st.write(df2)
            smile = df2['smiles'].values

            if st.sidebar.button("😊 GET STRUCTURE"):
                st.subheader('2-D structure')
                # raw_html = mols2grid.display(df2, mapping={"smiles": "SMILES"}, subset=["img", "iupac_name"], tooltip=["molecular_formula", "inchikey"])._repr_html_()
                raw_html = mols2grid.display(df2, mapping={"smiles": "SMILES"},
                                             tooltip=["SMILES"])._repr_html_()
                components.html(raw_html, width=900, height=900, scrolling=True)
