# vu-bioinformatics-projects

blastp -query queries/P00698.fasta -db db/db.fasta -outfmt 10
psiblast -query queries/P00698.fasta -db db/db.fasta -outfmt 10

python3 run_local_blast.py -ids SCOP_selections.txt -db ./db/db.fasta -outfile ./results/output.txt -outpng ./results/blast.png -q queries/ -vblast blast
python3 run_local_blast.py -ids SCOP_selections.txt -db ./db/db.fasta -outfile ./results/output.txt -outpng ./results/psiblast.png -q queries/ -vblast psiblast

e_value < 0.002
    voor blast: 
        648
    voor psiblast:
        671
