import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Draw

# Page Configuration
st.set_page_config(
    page_title="CADD Pocket Tool Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS: Blue, White & Yellow Scientific Theme
st.markdown("""
    <style>
    .stApp {
        background-color: #FFFFFF;
    }
    .main-header {
        font-size: 2.2rem;
        font-weight: 800;
        color: #0F294A;
        margin-bottom: 0.25rem;
    }
    .pro-badge {
        font-size: 0.95rem;
        font-weight: 700;
        background-color: #FEF08A;
        color: #854D0E;
        border: 1px solid #FACC15;
        padding: 3px 11px;
        border-radius: 12px;
        vertical-align: middle;
        margin-left: 6px;
    }
    .sub-header {
        font-size: 1.05rem;
        color: #1E40AF;
        margin-bottom: 1.2rem;
        font-weight: 500;
    }
    .section-banner {
        background: linear-gradient(90deg, #EFF6FF 0%, #FEFCE8 100%);
        border-left: 5px solid #1D4ED8;
        border-bottom: 1px solid #FDE047;
        padding: 10px 14px;
        border-radius: 6px;
        margin-bottom: 12px;
    }
    div[data-testid="stMetric"] {
        background-color: #F8FAFC;
        border: 1px solid #DBEAFE;
        border-top: 3px solid #2563EB;
        padding: 12px 14px;
        border-radius: 8px;
    }
    div[data-testid="stDownloadButton"] > button {
        background-color: #FACC15;
        color: #0F294A;
        font-weight: 700;
        border: 1px solid #EAB308;
    }
    div[data-testid="stDownloadButton"] > button:hover {
        background-color: #EAB308;
        color: #0F294A;
    }
    </style>
""", unsafe_allow_html=True)

# Header Section
st.markdown(
    '<div class="main-header">🧬 CADD Pocket Tool <span class="pro-badge">PRO v2.0</span></div>',
    unsafe_allow_html=True
)
st.markdown(
    '<div class="sub-header">In-Silico Molecular Property & Drug-Likeness Analyzer (Lipinski & Veber Rules)</div>',
    unsafe_allow_html=True
)

# Preset Famous Drugs Dictionary
PRESET_DRUGS = {
    "Custom SMILES Input": "",
    "Aspirin (Anti-inflammatory)": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine (Stimulant)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Paracetamol / Acetaminophen": "CC(=O)NC1=CC=C(O)C=C1",
    "Ibuprofen (Analgesic)": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Remdesivir (Antiviral - High MW)": "CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4"
}

# Input Section
col_input1, col_input2 = st.columns([1, 2])
with col_input1:
    selected_drug = st.selectbox("🧪 Quick Load Reference Compound:", list(PRESET_DRUGS.keys()))

with col_input2:
    default_smiles = PRESET_DRUGS[selected_drug] if selected_drug != "Custom SMILES Input" else "CC(=O)OC1=CC=CC=C1C(=O)O"
    smiles_input = st.text_input(
        "🔬 Enter Canonical / Isomeric SMILES String:",
        value=default_smiles,
        placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O"
    )

st.divider()

if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input.strip())
    
    if mol is not None:
        # Calculate Physicochemical Descriptors
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
        rot_bonds = Lipinski.NumRotatableBonds(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        
        # Lipinski's Rule of 5 Evaluation
        lipinski_violations = 0
        if mw > 500: lipinski_violations += 1
        if logp > 5: lipinski_violations += 1
        if hbd > 5: lipinski_violations += 1
        if hba > 10: lipinski_violations += 1
        
        # Veber's Rule Evaluation (Oral Bioavailability)
        veber_pass = (tpsa <= 140) and (rot_bonds <= 10)
        
        # Main Dashboard Layout (3 Columns)
        col1, col2, col3 = st.columns([1.1, 1.3, 1.3])
        
        # Column 1: 2D Molecular Structure
        with col1:
            st.markdown('<div class="section-banner"><strong>🖼️ 2D Structure</strong></div>', unsafe_allow_html=True)
            img = Draw.MolToImage(mol, size=(350, 350))
            st.image(img, use_container_width=True)
            st.caption(f"**Molecular Formula:** `{formula}`")
            
        # Column 2: Calculated Properties & Rules
        with col2:
            st.markdown('<div class="section-banner"><strong>📊 Physicochemical Profile</strong></div>', unsafe_allow_html=True)
            
            m1, m2 = st.columns(2)
            m1.metric("Molecular Weight", f"{mw} g/mol", delta="Pass (≤500)" if mw <= 500 else "High (>500)", delta_color="normal" if mw <= 500 else "inverse")
            m2.metric("LogP (Lipophilicity)", f"{logp}", delta="Pass (≤5)" if logp <= 5 else "High (>5)", delta_color="normal" if logp <= 5 else "inverse")
            
            m3, m4 = st.columns(2)
            m3.metric("H-Bond Donors", f"{hbd}", delta="Pass (≤5)" if hbd <= 5 else "High (>5)", delta_color="normal" if hbd <= 5 else "inverse")
            m4.metric("H-Bond Acceptors", f"{hba}", delta="Pass (≤10)" if hba <= 10 else "High (>10)", delta_color="normal" if hba <= 10 else "inverse")
            
            m5, m6 = st.columns(2)
            m5.metric("TPSA (Å²)", f"{tpsa} Å²", delta="Optimal (≤140)" if tpsa <= 140 else "High (>140)", delta_color="normal" if tpsa <= 140 else "inverse")
            m6.metric("Rotatable Bonds", f"{rot_bonds}", delta="Optimal (≤10)" if rot_bonds <= 10 else "Flexible (>10)", delta_color="normal" if rot_bonds <= 10 else "inverse")

        # Column 3: Radar Chart (Bioavailability Radar)
        with col3:
            st.markdown('<div class="section-banner"><strong>🕸️ Bioavailability Radar</strong></div>', unsafe_allow_html=True)
            
            # Normalize values relative to upper limits (1.0 = threshold limit)
            categories = ['MW (≤500)', 'LogP (≤5)', 'HBD (≤5)', 'HBA (≤10)', 'TPSA (≤140)', 'RotBonds (≤10)']
            normalized_vals = [
                min(mw / 500.0, 1.5),
                min(max(logp, 0) / 5.0, 1.5),
                min(hbd / 5.0, 1.5),
                min(hba / 10.0, 1.5),
                min(tpsa / 140.0, 1.5),
                min(rot_bonds / 10.0, 1.5)
            ]
            
            fig = go.Figure()
            # Threshold Boundary (Yellow/Gold Zone)
            fig.add_trace(go.Scatterpolar(
                r=[1, 1, 1, 1, 1, 1],
                theta=categories,
                fill='toself',
                fillcolor='rgba(250, 204, 21, 0.22)',
                line=dict(color='#EAB308', dash='dash', width=2),
                name='Optimal Drug Zone'
            ))
            # Compound Profile (Royal Blue Zone)
            fig.add_trace(go.Scatterpolar(
                r=normalized_vals,
                theta=categories,
                fill='toself',
                fillcolor='rgba(37, 99, 235, 0.28)',
                line=dict(color='#1D4ED8', width=2.5),
                name='Your Compound'
            ))
            fig.update_layout(
                paper_bgcolor='#FFFFFF',
                plot_bgcolor='#FFFFFF',
                polar=dict(
                    bgcolor='#F8FAFC',
                    radialaxis=dict(visible=True, range=[0, 1.5], showticklabels=False)
                ),
                showlegend=True,
                legend=dict(orientation="h", yanchor="bottom", y=-0.25, xanchor="center", x=0.5),
                margin=dict(l=30, r=30, t=20, b=40),
                height=310
            )
            st.plotly_chart(fig, use_container_width=True)

        st.divider()

        # Drug-Likeness Verdict & Report Export
        v_col1, v_col2, v_col3 = st.columns([1.5, 1.5, 1])
        
        with v_col1:
            if lipinski_violations <= 1:
                st.success(f"**✅ Lipinski's Rule of 5: PASSED** ({lipinski_violations} violation{'s' if lipinski_violations != 1 else ''})\n\nSuitable oral drug-likeness profile.")
            else:
                st.error(f"**❌ Lipinski's Rule of 5: FAILED** ({lipinski_violations} violations)\n\nPotential oral absorption/permeability challenges.")
                
        with v_col2:
            if veber_pass:
                st.success("**✅ Veber's Rule: PASSED**\n\nGood predicted oral bioavailability (TPSA ≤ 140 Å² & RotBonds ≤ 10).")
            else:
                st.warning("**⚠️ Veber's Rule: VIOLATED**\n\nHigh polar surface area or molecular flexibility may limit bioavailability.")

        with v_col3:
            # Prepare CSV Data for Download
            report_df = pd.DataFrame([{
                "SMILES": smiles_input,
                "Formula": formula,
                "Molecular_Weight": mw,
                "LogP": logp,
                "H_Bond_Donors": hbd,
                "H_Bond_Acceptors": hba,
                "TPSA": tpsa,
                "Rotatable_Bonds": rot_bonds,
                "Lipinski_Violations": lipinski_violations,
                "Lipinski_Status": "PASS" if lipinski_violations <= 1 else "FAIL",
                "Veber_Status": "PASS" if veber_pass else "FAIL"
            }])
            csv_data = report_df.to_csv(index=False).encode('utf-8')
            st.write("📥 **Export Research Data**")
            st.download_button(
                label="Download CSV Report",
                data=csv_data,
                file_name=f"CADD_Report_{formula}.csv",
                mime="text/csv",
                use_container_width=True
            )

    else:
        st.error("❌ Invalid SMILES string! Please check the chemical notation and try again.")

# Footer
st.markdown("---")
st.markdown("<div style='text-align: center; color: #1E3A8A; font-size: 0.85rem;'>Developed for Computational Drug Discovery & Cheminformatics Research | Powered by RDKit & Streamlit</div>", unsafe_allow_html=True)
