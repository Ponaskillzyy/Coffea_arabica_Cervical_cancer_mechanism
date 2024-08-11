# Mechanistic-Insights-into-the-Action-of-Coffea.-arabica-in-Cervical-Cancer-CC-Treatment.
In this study, we aimed to elucidate the mechanism of action of Coffea arabica (C. arabica) compounds in cervical cancer (CC) through a network-based approach. Our investigation commenced with a thorough examination of the bioactive compounds from C. arabica using the IMPPAT phytochemical screening database. Given the known issues of side effects and systemic toxicity associated with many small-molecule drugs, we prioritized evaluating the pharmacokinetic (PK) properties of these compounds. Our analysis led to the identification of seven compounds with favorable PK profiles, suggesting their potential for therapeutic application.

Subsequently, we utilized the Pharmapper database to map these compounds to their potential protein targets, creating a comprehensive list of C. arabica-targeted proteins. To contextualize these findings within the framework of cervical cancer, we analyzed bulk RNA sequencing data to identify differentially expressed genes (DEGs) in primary tumors (PT) relative to normal tissue (NT). This analysis revealed a total of 849 DEGs, providing a robust dataset for subsequent target identification.

Intersecting the DEGs with the list of compound-targeted proteins, we identified 38 core targets that overlap between C. arabica compounds and CC. This subset of genes represents potential targets modulated by C. arabica compounds in the context of cervical cancer.

To understand the biological significance of these core targets, we performed pathway enrichment analysis using the Enrichr R package and a curated list of 1,665 gene sets from the MSigDB Database, including Hallmark Gene Sets and Reactome canonical pathways. The enrichment analysis highlighted a significant association with 12 gene sets related to cell proliferation and survival. Notably, the core targets were highly enriched in gene sets associated with p53-mediated regulation of the cell cycle. This enrichment suggests that C. arabica compounds may influence CC cell proliferation and survival through modulation of the p53 pathway, a key regulator of the cell cycle and tumorigenesis.

In addition to the p53 pathway, our analysis also identified enrichment in other pathways known to be dysregulated in cancer, including SUMOylation, mTOR signaling, and MYC signaling. The integration of these results into a compound-target-pathway network revealed shared targets across these pathways, indicating the potential of C. arabica compounds to affect multiple oncogenic pathways simultaneously.

Given the strong association of the p53-mediated cell cycle regulation gene set with our core targets, we further investigated the prognostic implications of key genes within this set. We examined the correlation between the expression of these genes—TPX2, AURKA, CCNA2, MDM2, and CDK2—and patient survival outcomes in the TCGA dataset. Our analysis revealed that among these five genes, only CCNA2 was significantly associated with poor prognosis in cervical cancer patients. This finding positions CCNA2 as a potential biomarker for cervical cancer and suggests that it may be a critical target modulated by C. arabica compounds.
Cyclin A2, encoded by the CCNA2 gene, is a pivotal player in the regulation of the cell cycle, facilitating transitions through the G1/S and G2/M phases by activating cyclin-dependent kinase 2 (CDK2). The role of CCNA2 in cancer biology has been well-documented, particularly its contribution to uncontrolled cell proliferation and tumorigenesis. In the context of cervical cancer, CCNA2 is frequently upregulated, as shown in our study and other studies that correlate its high expression with enhanced tumor growth and progression. This upregulation is thought to promote cell cycle progression, thereby enabling the rapid proliferation of cancer cells and contributing to the aggressive nature of the disease.
Given the overexpression of CCNA2 in cervical cancer, it holds significant promise as a prognostic biomarker. Its elevated levels have been associated with advanced tumor stages and poor clinical outcomes, underscoring its potential utility in predicting disease prognosis. Moreover, targeting CCNA2 may offer therapeutic benefits, as inhibiting its activity could disrupt the aberrant cell cycle progression characteristic of cancer cells.
We hypothesize that bioactive compounds derived from C. arabica may potentially inhibit CCNA2 in cervical cancer cells, leading to cell cycle arrest, reduced proliferation, and enhanced apoptosis, ultimately limiting tumor growth. To explore this hypothesis, we performed molecular docking studies to assess the binding affinities of seven compounds from Coffea arabica against CCNA2, using a standard cell cycle inhibitor (Tagtociclib) known to promote apoptosis as a control. Our study revealed that stigmasterol had better predicted binding affinity than the control drug. Additionally, we computationally predict the pIC50 values of stigmasterol and the control using a machine learning model trained on cell proliferation inhibitors with known IC50 values. We observed similar predicted pIC50 values between stigmasterol and the control, suggesting its potential as an inhibitor of cervical cancer cell proliferation. 
In summary, our study suggests that compounds from Coffea arabica, particularly stigmasterol, may inhibit CCNA2 in cervical cancer cells, thereby modulating the p53-mediated cell cycle pathway to arrest cell proliferation and induce apoptosis.



