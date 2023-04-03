from flask import Flask, request, render_template
from pep import PeptideOperations
import pep
from Bio import SeqIO
import io



app = Flask(__name__, template_folder='./template', static_folder='./static')
app.config["DEBUG"] = True


@app.route("/")
def home():
    return render_template('home.html')


@app.route("/peptool")
def peptool():
    return render_template('pepTool.html')


@app.route("/blast_tool")
def blasttool():
    return render_template('blastTool.html')


@app.route("/calculate_result", methods=['POST'])
def calculate():
    pepOps = PeptideOperations()
    seq = None
    check_list= []
    result=None

    seq = request.form["seq"]
    check_list=request.form.getlist("cal")

    result = pepOps.check(seq, check_list)
    return render_template('pepResult.html', result=result )


@app.route("/calculate_blast_result", methods=['POST'])
def calculateBlast():
    result=None
    fileContent=None
    # Read the uploaded file
    fasta_file = request.files.get('fasta_file')
    if fasta_file is not None and fasta_file.filename != '':
        contents = fasta_file.read()
        # Convert the contents to a string
        fileContent = SeqIO.read(io.StringIO(contents.decode('utf-8')), 'fasta')
        result = pep.blast(fileContent)
    return render_template('blastResult.html', result=result)


@app.route("/documentation")
def documentation():
    return render_template('help.html')


if __name__ == '__main__':
    app.run()