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


# Assignment GO

## Notes

Go-terms for a single protein can be retrieved by accessing `http://www.uniprot.org/uniprot/?query=accession:<PROTID>&format=tab&columns=id,organism-id,go-id,pathway` and by replacing `<PROTID>` with the protein ID.

## Question 3

Convergent evolution...

## Question 6

### PART A

Go-terms for protein **LAP3** - P28838 (AMPL_HUMAN):
```
Entry	Organism ID	Gene ontology IDs	Pathway
Q01532	559292	GO:0000122; GO:0000978; GO:0003690; GO:0003697; GO:0003729; GO:0004197; GO:0005737; GO:0005739; GO:0008234; GO:0009636; GO:0043418; GO:0046677	
```

Go-terms for protein **BLMH** - Q13867 (BLMH_HUMAN):
```
Entry	Organism ID	Gene ontology IDs	Pathway
Q13867	9606	GO:0000209; GO:0004177; GO:0004180; GO:0004197; GO:0005634; GO:0005737; GO:0005829; GO:0006508; GO:0008234; GO:0009636; GO:0042493; GO:0042802; GO:0043418; GO:0070062	
```

lap3 $\cap$ blmh:
```python
>>> lap3 = set("GO:0000122; GO:0000978; GO:0003690; GO:0003697; GO:0003729; GO:0004197; GO:0005737; GO:0005739; GO:0008234; GO:0009636; GO:0043418; GO:0046677".split(';')) 
>>> blmh = set("GO:0000209; GO:0004177; GO:0004180; GO:0004197; GO:0005634; GO:0005737; GO:0005829; GO:0006508; GO:0008234; GO:0009636; GO:0042493; GO:0042802; GO:0043418; GO:0070062".split(';'))

# intersection = g(p1) intersection g(p2)
>>> intersection = lap3.intersection(blmh)
{' GO:0004197', ' GO:0005737', ' GO:0008234', ' GO:0009636', ' GO:0043418'}
# len(intersection)
>>> len(intersection)
5

# union = g(p1) union g(p2)
>>> union = lap3.union(blmh)
{' GO:0005737', ' GO:0005634', ' GO:0070062', ' GO:0004180', ' GO:0005829', ' GO:0009636', ' GO:0005739', ' GO:0004197', ' GO:0006508', ' GO:0003690', ' GO:0003697', 'GO:0000209', ' GO:0003729', ' GO:0042802', ' GO:0004177', ' GO:0000978', 'GO:0000122', ' GO:0046677', ' GO:0043418', ' GO:0008234', ' GO:0042493'}
>>> len(union)
21

# |intersection| / |union|
>>> float(len(intersection))/len(union)
0.23809523809523808
```

### PART B

The score of a protein s(p1, p2) where p1==p2 would be 1. Proof:

Within set theory, the idempotency law states that any union or intersection of any set with itself is equal to itself.
$A \cap A = A$
$A \cup A = A$

therefore:
$$s(A, A)=\frac{\left |A \cap A  \right |}{\left |A \cup A  \right |} = \frac{\left | A \right |}{\left | A \right |} = 1$$

## Question 7

To determine the number of unique combinations between two equal sets of proteins, the following formula applies.(n*(n-1))/2The total is divided by two because identical pairs would be included assuming that any combination (as opposed to a permutation) is seen as equal (eg: (a,b)=(b,a)).In the case of n being 200, the total number of unique combinations with same-protein pairs being excluded is:(200*199)/2 = 19900So. a total of 19990 pairs would be scored for.

# Q

1. Combinaties of permutaties?
2. bereik similarity inclusive of **exclusive** range