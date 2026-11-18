import streamlit as st

# Page Configuration
st.set_page_config(
    page_title="CADD Pocket Tool",
    page_icon="💊",
    layout="centered"
)

# App Title & Description
st.title("💊 CADD Pocket Tool")
st.subheader("Lipinski's Rule of 5 Analyzer for Researchers")
st.write("Enter the molecular properties below to evaluate the drug-likeness of your compound.")

st.markdown("---")

# User Inputs using Streamlit widgets (Great for mobile & desktop)
mw = st.number_input("Molecular Weight (g/mol)", min_value=0.0, max_value=2000.0, value=180.16, step=0.1)
logp = st.number_input("LogP (Lipophilicity)", min_value=-10.0, max_value=20.0, value=1.2, step=0.1)
donors = st.number_input("H-Bond Donors", min_value=0, max_value=50, value=1, step=1)
acceptors = st.number_input("H-Bond Acceptors", min_value=0, max_value=50, value=3, step=1)

st.markdown("---")

# Evaluate Button
if st.button("Evaluate Drug-Likeness", type="primary"):
    violations = 0
    results = []
    
    # 1. Molecular Weight Check
    if mw <= 500:
        results.append(("Molecular Weight", mw, "Passed (<= 500 g/mol)", True))
    else:
        results.append(("Molecular Weight", mw, "Violation (> 500 g/mol)", False))
        violations += 1
        
    # 2. LogP Check
    if logp <= 5:
        results.append(("LogP", logp, "Passed (<= 5)", True))
    else:
        results.append(("LogP", logp, "Violation (> 5)", False))
        violations += 1
        
    # 3. H-Bond Donors Check
    if donors <= 5:
        results.append(("H-Bond Donors", donors, "Passed (<= 5)", True))
    else:
        results.append(("H-Bond Donors", donors, "Violation (> 5)", False))
        violations += 1
        
    # 4. H-Bond Acceptors Check
    if acceptors <= 10:
        results.append(("H-Bond Acceptors", acceptors, "Passed (<= 10)", True))
    else:
        results.append(("H-Bond Acceptors", acceptors, "Violation (> 10)", False))
        violations += 1

    # Display Results in UI
    st.markdown("### 📊 Analysis Report")
    
    for param, val, desc, passed in results:
        if passed:
            st.success(f"**{param}** ({val}): {desc}")
        else:
            st.error(f"**{param}** ({val}): {desc}")
            
    st.info(f"**Total Violations:** {violations}")
    
    # Final Verdict
    if violations == 0:
        st.success("🎉 **RESULT:** Excellent Drug-Likeness (0 Violations). Highly viable candidate!")
    elif violations == 1:
        st.warning("⚠️ **RESULT:** Good Drug-Likeness (1 Violation permitted under extended rules).")
    else:
        st.error("❌ **RESULT:** Poor Drug-Likeness (Multiple violations. High failure risk in oral bioavailability).")