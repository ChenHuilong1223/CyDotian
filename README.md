# **CyDotian**’s manual

A toolkit for identifying internal repeats of nucleic acid and amino acid sequences.

| Authors | Huilong Chen           |
| ------- | ---------------------- |
| Email   | chenhuilong131@163.com |

# Overview

CyDotian algorithm is a dynamic programming algorithm, which can identify all internal repeats of the sequence itself that allow mutation. In order to achieve efficient output of the results, we used C, the fastest underlying computer language available, to implement the algorithm and allow it to be compiled into an executable program. Downstream analysis tools are then written in Python, the most popular language for data processing and visualisation. These downstream analysis tools include processing the location and number of repeat segments, plotting dotplots, plotting depth plots, calculating repetition density and outputting specific repeat segment comparison details. All are batch processed and exported, which is extremely user-friendly. Users can suggest and optimise the development of all codes. Moreover, due to the applicability of the CyDotian algorithm, it can also be used to identify similarity segments between two different sequences that allow for mutations. In addition, by replacing the U in the RNA sequence with a T, CyDotian can be used to find all the reverse complementary sequences, helping to predict the stem-loop structure of the RNA molecule.



Usage information is built into the program. To display usage on the screen, the user simply runs the program by specifying the -h/--help parameter:

$ python3 program_name.py -h/--help (for Python scripts)



###### The following is the list of executable programs:

bpRepeatScan (used for nucleic acid sequences)

aaRepeatScan (used for amino acid sequences)

slidingWindow (used for used for nucleic acid and amino acid sequences via sliding window method)



###### Parameter configuration file

CyDotian.config (used for CyDotian algorithm)

CyDotian_exact_match.config (used for MUMmer's repeat-match algorithm)

CyDotian_sliding_window.config (used for sliding window algorithm)



###### Batch processing programs:

Tool 1. 1.0_batch_check_sequence.py

Tool 2. 1.1_batch_run_CyDotian.py

Tool 3. 1.2_batch_run_draw_dotplot.py

Tool 4. 1.3_batch_run_draw_depth_plot.py

Tool 5. 1.4_batch_run_output_repeat_density.py

Tool 6. 1.5_batch_extract_repeat_sequences.py

Tool 7. 1.6_Extract_the_corresponding_results_by_name.py

Tool 8. 1.7_batch_run_CyDotian_in_pairwise_comparison_mode.py

Tool 9. 1.8_batch_extract_repeat_sequences_in_pairwise_comparison_mode.py

Tool 10. 2.1_batch_run_CyDotian_exact_match.py

Tool 11. 2.7_batch_run_CyDotian_exact_match_in_pairwise_comparison_mode.py

Tool 12. 3.1_batch_run_CyDotian_sliding_window.py

Tool 13. 3.2_batch_run_CyDotian_sliding_window_in_pairwise_comparison_mode.py

Tool 14. 3.3_batch_run_draw_dotplot_sliding_window.py



###### Articles that have cited the CyDotian algorithm:

*Ge W, Chen H, Zhang Y, et al. Integrative genomics analysis of the ever-shrinking pectin methylesterase (PME) gene family in foxtail millet (Setaria italica)[J]. Functional Plant Biology, 2022, 49(10): 874-886.. https://doi.org/10.1071/FP21319*

# Documentation

The online documentation is located at the [GitHub Wiki](https://github.com/ChenHuilong1223/CyDotian/wiki).
