# susD screen

Check of the annotated fungal proteomes for homologs of the *Agrobacterium* succinamopine-synthase-associated genes
susD (ARU12249.1) and susL (AAZ50401.1) of pTiEU6; both queries are in `susD.faa`.

* `screen_one.sh`: DIAMOND blastp (`--very-sensitive`, bit score >= 50), one proteome per call, over the same proteome list as the main screen
  (`analysis/refset/screen/step3/proteomes.txt`). Run it with LSF job arrays (`LSB_JOBINDEX`) or directly (`screen_one.sh 0 <list>`).
* `best_hits.tsv`: best hit per proteome (seed, assembly, organism, protein, identity, bit score, alignment length, seed length).
* `SUMMARY.txt`: result in two lines.
* The per-proteome output folder `out/` (one TSV per proteome) is not copied into the repository; it is recreated by `screen_one.sh`.
