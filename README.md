# PepAnalyzer
PepAnalyzer is a powerful and user-friendly tool designed to assist researchers and students in analysing and understanding various aspects of peptides. With this tool, you can easily calculate molecular weight, predict secondary structure, predict hydrophobicity, and much more, making it a helpful tool for students and researchers to get insight about peptides.




## Features
Input a peptide sequence, either in upper or lowercase letters. The app will show an alert if ambiguity is detected in the sequence.
There are 16 types of functionalities available to be calculated for the given peptide sequence and an extra added Blast tool access, including:
1. Physical/Chemical Properties:
  - Peptide Size
  - Molecular Weight
  - Distribution of Amino Acid
  - Amino Acid Classifier
2. Predictive Tools:
 - Peptide Charge Calculator
 - Aromicity Calculator
 - Half Life Predictor
 - Predict Secondary Structure
 - Molar Extinction Coefficient
 - Disulfide Bridges Calculation
 - Binding Potential Calculator
3. Graphical Representations:
 - Hydropathy Plot for Peptide
 - Helical Wheel Projection for Peptide
 - Distribution of Amino Acid
 - Amino Acid Classifier
4. Conversion Utility:
 - One-to-Three Letter Code Conversion
5. Bio Blast tool




## Installation

Run the following code in terminal to install necessery libraries.
```html
pip install -r requirements.txt
```
modlamp==4.3.0
biopython==1.79
mpld3==0.5.5
numpy==1.21.5
Flask==2.1.1




## Running the webapp

Run the following code in terminal to use the webapp locally.
```html
python app.py
```




## Help and Documentation
To learn more about how to use PepAnalyzer, visit the Help and Documentation section in the application. Here's an overview of what you can expect to find:

Sequence: Input a peptide sequence, either in upper or lowercase letters. The app will show an alert if ambiguity is detected in the sequence.
Set parameters: There are 16 types of functionalities available to be calculated for the given peptide sequence and an extra added Blast tool access.
Submit: Clicking on Submit will take us to a new page where we can see all the results that have been calculated. If the blast tool was accessed, this may take some time to display results depending on the length of the input sequence.

Query: Mail your query to [mrlabdbc@db.du.ac.in](mailto:mrlabdbc@db.du.ac.in) about paper.

Query: Mail your query to [aswanthmanoj51@gmail.com](mailto:aswanthmanoj51@gmail.com) about code.
