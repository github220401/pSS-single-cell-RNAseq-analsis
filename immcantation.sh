#immcantation
#https://github.com/Meng0625/immcantation
source activate changeo

# Linux/Mac OS X
docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.3.0 bash

# Run 10x Genomics processing script
changeo-10x -s filtered_contig.fasta -a filtered_contig_annotations.csv -o . \
        -g human -t ig -x 0.1

MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s filtered_contig.fasta \
                   -r /home/database/human_germline/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/*.fasta --10x filtered_contig_annotations.csv --extended

ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IGH" \
                           --logic all --regex --outname heavy
ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IG[LK]" \
                           --logic all --regex --outname light

# define heavy chain clones
DefineClones.py -d heavy_parse-select.tsv --act set --model ham \
       --norm len --dist 0.09 --outname filtered_contig_heavy

light_cluster.py -d filtered_contig_heavy_clone-pass.tsv -e light_parse-select.tsv \
               -o 10X_s0_clone-pass.tsv
#light_cluster.py -d filtered_contig_heavy_clone-pass.tsv -e light_parse-select.tsv \
#                              -o 10X_s1_clone-pass.tsv
#light_cluster.py -d filtered_contig_heavy_clone-pass.tsv -e light_parse-select.tsv \
#                              -o 10X_s2_clone-pass.tsv
#light_cluster.py -d filtered_contig_heavy_clone-pass.tsv -e light_parse-select.tsv \
#                              -o 10X_hc_clone-pass.tsv

# split heavy chain clones with different light chains
light_cluster.py -d filtered_contig_heavy_clone-pass.tsv -e filtered_contig_light_productive-T.tsv \
    -o filtered_contig_heavy_clone-light.tsv

# reconstruct heavy chain germline V and J sequences
CreateGermlines.py -d filtered_contig_heavy_clone-light.tsv -g dmask --cloned \
        -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta \
        /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta \
        /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta \
        --outname filtered_contig_heavy

BuildTrees.py -d filtered_contig_heavy_germ-pass.tsv --minseq 3 --clean all \
    --igphyml --collapse --nproc 2 --asr 0.1
