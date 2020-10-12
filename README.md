# vu-bioinformatics-projects

blastp -query queries/P00698.fasta -db db/db.fasta -outfmt 10
psiblast -query queries/P00698.fasta -db db/db.fasta -outfmt 10

python3 run_local_blast.py -ids SCOP_selections.txt -db ./db/db.fasta -outfile ./results/output_blast.txt -outpng ./results/blast.png -q queries/ -vblast blast
python3 run_local_blast.py -ids SCOP_selections.txt -db ./db/db.fasta -outfile ./results/output_psiblast.txt -outpng ./results/psiblast.png -q queries/ -vblast psiblast

e_value < 0.002
    voor blast: 
        648
    voor psiblast:
        671

python3 classify_go.py -uniprot SCOP_selections.txt -threshold1 0.054 -threshold2 0.111 -output_file results/go_output.txt

python3 roc_plot.py -blast_results ./results/output.txt -go_results ./results/go_output.txt -outpng roc.png

