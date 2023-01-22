#!/usr/bin/env python3
# John Vidale 4/2020
# last modified 8/2022

def run_compare_ind(repeater = 'NoName',do_global = False, do_YKA = False, do_ILAR = False):

    import os
    import sys
    import time
    import matplotlib.pyplot as plt
    from termcolor import colored

    #%% Import functions
    pro_directory = '/Users/vidale/Documents/GitHub/Array_codes/Run_YKA_WRA'
    os.chdir(pro_directory)
    from run_compare_global     import run_compare_global
    from run_individual_df      import run_individual_df
    from run_compare_pair       import run_compare_pair

    start_time_wc = time.time()

    plt.close('all')

    if   repeater == 'P01':
        eq_num1 = 701; eq_num2 = 726; do_global = False; do_ILAR = False; do_YKA =  True; tshift = -5.77; shift_both =  3.40; flip =  True
    elif repeater == 'P02':
        eq_num1 = 701; eq_num2 = 753; do_global = False; do_ILAR = False; do_YKA =  True; tshift = -5.00; shift_both =  3.40; flip =  True
    elif repeater == 'P03':
        eq_num1 = 702; eq_num2 = 708; do_global = False; do_ILAR = False; do_YKA =  True; tshift = -7.98; shift_both =  2.00; flip =  True
    elif repeater == 'P04':
        eq_num1 = 702; eq_num2 = 723; do_global = False; do_ILAR = False; do_YKA =  True; tshift = -6.70; shift_both =  2.50; flip =  True
    elif repeater == 'P05':
        eq_num1 = 702; eq_num2 = 739; do_global = False; do_ILAR = False; do_YKA =  True; tshift = -4.78; shift_both =  3.00; flip =  True
    elif repeater == 'P06':
        eq_num1 = 702; eq_num2 = 757; do_global = False; do_ILAR = False; do_YKA =  True; tshift = -4.91; shift_both =  3.00; flip =  True
    elif repeater == 'P07':
        eq_num1 = 703; eq_num2 = 721; do_global =  True; do_ILAR = False; do_YKA =  True; tshift =  0.00; shift_both = -1.30; flip =  True
    elif repeater == 'P08':
        eq_num1 = 704; eq_num2 = 706; do_global =  True; do_ILAR = False; do_YKA =  True; tshift = -5.35; shift_both =  7.00; flip = False
    elif repeater == 'P09':
        eq_num1 = 704; eq_num2 = 713; do_global =  True; do_ILAR = False; do_YKA =  True; tshift = -3.43; shift_both =  5.20; flip = False
    elif repeater == 'P10':
        eq_num1 = 705; eq_num2 = 711; do_global =  True; do_ILAR = False; do_YKA =  True; tshift =  0.04; shift_both =  0.00; flip = False
    elif repeater == 'P11':
        eq_num1 = 705; eq_num2 = 729; do_global =  True; do_ILAR = False; do_YKA =  True; tshift = -3.85; shift_both = -0.50; flip = False
    elif repeater == 'P12':
        eq_num1 = 706; eq_num2 = 713; do_global =  True; do_ILAR = False; do_YKA =  True; tshift =  1.95; shift_both =  0.00; flip = False
    elif repeater == 'P13':
        eq_num1 = 707; eq_num2 = 756; do_global =  True; do_ILAR = False; do_YKA =  True; tshift =  1.80; shift_both =  0.00; flip = False
    elif repeater == 'P14':
        eq_num1 = 709; eq_num2 = 730; do_global =  True; do_ILAR = False; do_YKA =  True; tshift =  2.28; shift_both = -2.00; flip = False
    elif repeater == 'P15':
        eq_num1 = 709; eq_num2 = 744; do_global =  True; do_ILAR = False; do_YKA =  True; tshift = -2.59; shift_both =  0.00; flip = False
    elif repeater == 'P16':
        eq_num1 = 708; eq_num2 = 723; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  1.29; shift_both = -5.00; flip =  True
    elif repeater == 'P17':
        eq_num1 = 708; eq_num2 = 739; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  3.22; shift_both = -4.50; flip =  True
    elif repeater == 'P18':
        eq_num1 = 708; eq_num2 = 757; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  3.07; shift_both = -6.00; flip =  True
    elif repeater == 'P19':
        eq_num1 = 710; eq_num2 = 728; do_global =  True; do_ILAR = False; do_YKA =  True; tshift = -0.60; shift_both =  0.00; flip = False
    elif repeater == 'P20':
        eq_num1 = 711; eq_num2 = 729; do_global =  True; do_ILAR = False; do_YKA =  True; tshift = -3.89; shift_both =  0.00; flip = False
    elif repeater == 'P21':
        eq_num1 = 712; eq_num2 = 718; do_global =  True; do_ILAR = False; do_YKA =  True; tshift =  0.35; shift_both =  0.00; flip = False
    elif repeater == 'P22':
        eq_num1 = 714; eq_num2 = 741; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -2.90; shift_both =  1.80; flip = False
    elif repeater == 'P23':
        eq_num1 = 714; eq_num2 = 748; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  1.95; shift_both =  1.80; flip = False
    elif repeater == 'P24':
        eq_num1 = 715; eq_num2 = 758; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -3.60; shift_both = -1.18; flip =  True
    elif repeater == 'P25':
        eq_num1 = 716; eq_num2 = 737; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -2.24; shift_both =  0.00; flip =  True
    elif repeater == 'P26':
        eq_num1 = 717; eq_num2 = 743; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -4.98; shift_both =  0.28; flip = False
    elif repeater == 'P27':
        eq_num1 = 719; eq_num2 = 754; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  0.03; shift_both = -0.10; flip =  True
    elif repeater == 'P28':
        eq_num1 = 720; eq_num2 = 736; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -3.95; shift_both =  0.05; flip = False
    elif repeater == 'P29':
        eq_num1 = 722; eq_num2 = 740; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -2.30; shift_both = -0.12; flip =  True
    elif repeater == 'P30':
        eq_num1 = 723; eq_num2 = 739; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  1.94; shift_both = -3.72; flip =  True
    elif repeater == 'P31':
        eq_num1 = 723; eq_num2 = 757; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  1.78; shift_both = -3.72; flip =  True
    elif repeater == 'P32':
        eq_num1 = 724; eq_num2 = 746; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  0.75; shift_both = -1.14; flip = False
    elif repeater == 'P33':   # weird that global disagrees, but ILAR is very similar
        eq_num1 = 725; eq_num2 = 749; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  0.56; shift_both = -2.15; flip = False
    elif repeater == 'P34':
        eq_num1 = 726; eq_num2 = 753; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  0.74; shift_both = -2.09; flip =  True
    elif repeater == 'P35':
        eq_num1 = 727; eq_num2 = 752; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  3.00; shift_both = -2.36; flip =  True
    elif repeater == 'P36':
        eq_num1 = 730; eq_num2 = 744; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -4.82; shift_both =  0.00; flip =  True
    elif repeater == 'P37':
        eq_num1 = 731; eq_num2 = 747; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -1.41; shift_both = -2.41; flip =  True
    elif repeater == 'P38':  # big changes only in YKA
        eq_num1 = 732; eq_num2 = 750; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  0.60; shift_both = -2.27; flip =  True
    elif repeater == 'P39':  # some surprising global shifts
        eq_num1 = 733; eq_num2 = 745; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -2.12; shift_both = -1.62; flip =  True
    elif repeater == 'P40':
        eq_num1 = 734; eq_num2 = 742; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  7.70; shift_both = -6.18; flip = False
    elif repeater == 'P41':
        eq_num1 = 735; eq_num2 = 751; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  2.58; shift_both = -1.85; flip = False
    elif repeater == 'P42':
        eq_num1 = 738; eq_num2 = 755; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -0.15; shift_both = -2.56; flip =  True
    elif repeater == 'P43':  # big changes only in YKA
        eq_num1 = 739; eq_num2 = 757; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -0.15; shift_both = -1.77; flip =  True
    elif repeater == 'P44':  # big changes only in YKA
        eq_num1 = 741; eq_num2 = 748; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  4.86; shift_both = -1.61; flip = False
    elif repeater == 'P45':  # Sumatra, IV2 is too close
        eq_num1 = 759; eq_num2 = 760; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift = -0.67; shift_both = -0.45; flip =  True
# 2nd set of repeaters
    elif repeater == 'P50':
        eq_num1 = 712; eq_num2 = 848; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -4.15; shift_both = 00; flip =  False
    elif repeater == 'P51':
        eq_num1 = 718; eq_num2 = 848; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -4.50; shift_both = 00; flip =  False
    elif repeater == 'P52':
        eq_num1 = 734; eq_num2 = 848; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   1.10; shift_both = -6.00; flip =  False
    elif repeater == 'P53':
        eq_num1 = 742; eq_num2 = 848; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -6.60; shift_both =  1.50; flip =  False
    elif repeater == 'P54':
        eq_num1 = 704; eq_num2 = 827; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift = -11.28; shift_both =  6.00; flip =  False
    elif repeater == 'P55':
        eq_num1 = 706; eq_num2 = 827; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -5.90; shift_both = 00; flip =  False
    elif repeater == 'P56':
        eq_num1 = 713; eq_num2 = 827; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -7.853; shift_both =  4.00; flip =  False
    elif repeater == 'P57':
        eq_num1 = 704; eq_num2 = 849; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -9.82; shift_both =  6.00; flip =  False
    elif repeater == 'P58':
        eq_num1 = 706; eq_num2 = 849; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -4.45; shift_both = 00; flip =  False
    elif repeater == 'P59':
        eq_num1 = 713; eq_num2 = 849; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -6.40; shift_both =  4.00; flip =  False
    elif repeater == 'P60':
        eq_num1 = 719; eq_num2 = 811; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   0.00; shift_both = 00; flip =  False
    elif repeater == 'P61':
        eq_num1 = 719; eq_num2 = 822; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   0.67; shift_both = 00; flip =  False
    elif repeater == 'P62':
        eq_num1 = 822; eq_num2 = 754; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -0.65; shift_both = 00; flip =  False
    elif repeater == 'P63':
        eq_num1 = 811; eq_num2 = 754; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   0.04; shift_both = 00; flip =  False
    elif repeater == 'P64':
        eq_num1 = 802; eq_num2 = 715; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -0.55; shift_both = 00; flip =  False
    elif repeater == 'P65':
        eq_num1 = 802; eq_num2 = 758; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -4.12; shift_both = 00; flip =  False
    elif repeater == 'P66':
        eq_num1 = 715; eq_num2 = 830; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -4.93; shift_both = 00; flip =  False
    elif repeater == 'P67':
        eq_num1 = 830; eq_num2 = 758; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   1.31; shift_both = 00; flip =  False
    elif repeater == 'P68':
        eq_num1 = 807; eq_num2 = 810; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -7.63; shift_both = -2.00; flip =  False
    elif repeater == 'P69':
        eq_num1 = 807; eq_num2 = 825; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -6.00; shift_both = 00; flip =  False
    elif repeater == 'P70':
        eq_num1 = 807; eq_num2 = 833; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -4.66; shift_both = 00; flip =  False
    elif repeater == 'P71':
        eq_num1 = 810; eq_num2 = 825; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   1.63; shift_both = -6.00; flip =  False
    elif repeater == 'P72':
        eq_num1 = 810; eq_num2 = 833; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   2.97; shift_both = -10; flip =  False
    elif repeater == 'P73':
        eq_num1 = 825; eq_num2 = 833; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   1.35; shift_both = -5.00; flip =  False
    elif repeater == 'P74':
        eq_num1 = 812; eq_num2 = 823; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   3.66; shift_both = 00; flip =  False
    elif repeater == 'P75':
        eq_num1 = 812; eq_num2 = 842; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   3.23; shift_both = -2.50; flip =  False
    elif repeater == 'P76':
        eq_num1 = 823; eq_num2 = 842; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -0.42; shift_both =  2.0; flip =  False
    elif repeater == 'P77':
        eq_num1 = 707; eq_num2 = 809; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =   1.60; shift_both = 00; flip =  False
    elif repeater == 'P78':
        eq_num1 = 809; eq_num2 = 756; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   0.04; shift_both =  1.50; flip =  False
    elif repeater == 'P79':
        eq_num1 = 720; eq_num2 = 845; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -2.98; shift_both = 00; flip =  False
    elif repeater == 'P80':
        eq_num1 = 736; eq_num2 = 845; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   0.99; shift_both = -3.9; flip =  False
    elif repeater == 'P81':
        eq_num1 = 808; eq_num2 = 824; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -6.49; shift_both = 00; flip =  False
    elif repeater == 'P82':
        eq_num1 = 808; eq_num2 = 837; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -6.98; shift_both = 00; flip =  False
    elif repeater == 'P83':
        eq_num1 = 824; eq_num2 = 837; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -0.48; shift_both = -6.00; flip =  False
    elif repeater == 'P84':
        eq_num1 = 801; eq_num2 = 727; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -2.20; shift_both = 00; flip =  False
    elif repeater == 'P85':
        eq_num1 = 801; eq_num2 = 752; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =   0.80; shift_both = 00; flip =  False
    elif repeater == 'P86':
        eq_num1 = 804; eq_num2 = 847; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -1.75; shift_both = 00; flip =  False
    elif repeater == 'P87':
        eq_num1 = 820; eq_num2 = 835; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -0.22; shift_both = 00; flip =  False
    elif repeater == 'P88':
        eq_num1 = 819; eq_num2 = 840; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   2.88; shift_both = 00; flip =  False
    elif repeater == 'P89':
        eq_num1 = 817; eq_num2 = 834; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   4.97; shift_both = -2.50; flip =  False
    elif repeater == 'P90':
        eq_num1 = 815; eq_num2 = 831; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   3.00; shift_both = 00; flip =  False
    elif repeater == 'P91':
        eq_num1 = 814; eq_num2 = 841; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   2.67; shift_both = 00; flip =  False
    elif repeater == 'P92':
        eq_num1 = 813; eq_num2 = 844; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   0.24; shift_both = 00; flip =  False
    elif repeater == 'P93':
        eq_num1 = 800; eq_num2 = 828; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -6.20; shift_both = 00; flip =  False
    elif repeater == 'P94':
        eq_num1 = 805; eq_num2 = 816; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -0.50; shift_both = 00; flip =  False
    elif repeater == 'P95':
        eq_num1 = 803; eq_num2 = 836; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =   1.92; shift_both = 00; flip =  False
    elif repeater == 'P96':
        eq_num1 = 714; eq_num2 = 832; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   1.02; shift_both = 00; flip =  False
    elif repeater == 'P97':
        eq_num1 = 806; eq_num2 = 818; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -9.43; shift_both = 00; flip =  False
    elif repeater == 'P98':
        eq_num1 = 821; eq_num2 = 838; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -1.00; shift_both = 00; flip =  False
    elif repeater == 'P99':
        eq_num1 = 826; eq_num2 = 846; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -0.36; shift_both = 00; flip =  False
    elif repeater == 'P100':
        eq_num1 = 829; eq_num2 = 843; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   3.04; shift_both = -2.50; flip =  False
    elif repeater == 'P101':
        eq_num1 = 839; eq_num2 = 850; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -0.84; shift_both = 00; flip =  False
    elif repeater == 'P102':
        eq_num1 = 712; eq_num2 = 742; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =   2.44; shift_both = 00; flip =  False
    elif repeater == 'P103':
        eq_num1 = 712; eq_num2 = 734; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -5.20; shift_both = 00; flip =  False
    elif repeater == 'P104':
        eq_num1 = 718; eq_num2 = 734; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =  -5.60; shift_both = 00; flip =  False
    elif repeater == 'P105':
        eq_num1 = 718; eq_num2 = 742; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   2.10; shift_both = 00; flip =  False
    elif repeater == 'P106':
        eq_num1 = 827; eq_num2 = 849; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   1.47; shift_both = -4.00; flip =  False
    elif repeater == 'P107':
        eq_num1 = 811; eq_num2 = 822; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   0.70; shift_both = 00; flip =  False
    elif repeater == 'P108':
        eq_num1 = 802; eq_num2 = 830; do_global =  True; do_ILAR =  False; do_YKA =  True; tshift =  -5.44; shift_both = 00; flip =  False
    elif repeater == 'P109':
        eq_num1 = 832; eq_num2 = 748; do_global =  True; do_ILAR =  True; do_YKA =  True; tshift =   1.00; shift_both = 00; flip =  False
    else:
        print(colored('No such event pair ' + repeater, 'yellow'))
        sys.exit(-1)

    # freq_min = 0.6; freq_max = 3
    freq_min = 1; freq_max = 2
    # freq_min = 2; freq_max = 4
    do_ILAR_pre = do_ILAR

    do_YKA = False
    do_ILAR = False
    # do_ILAR_pre = False
    do_global = False

    if do_global:
        start_buff = -20 # analysis window start relative to phase arrival
        wind_len    = 50 # analysis window length
        run_compare_global(repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2,
                freq_min = freq_min, freq_max = freq_max, ARRAY = 7, tshift = tshift,
                start_buff = start_buff, wind_len = wind_len)

    wind_buff = 10 # buffer before and after time window of analysis
    plot_peak = 1

    beam_offset = 0.02
    beam_width  = 0.015
    # slow_delta  = 0.004
    slow_delta  = 0.005

    if do_YKA:
        Zstart_buff = -40 # analysis window start relative to phase arrival
        wind_len    = 120 # analysis window length
        run_compare_pair(repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2, dphase = 'PKIKP', shift_both = shift_both,
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 5,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff, flip = flip,
                tshift = tshift, fig_index = 200, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

    if do_ILAR:
        Zstart_buff = -2 # analysis window start relative to phase arrival
        wind_len    =  7 # analysis window length
        run_compare_pair(repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2, dphase = 'PKIKP', shift_both = shift_both,
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 6,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff, flip = flip,
                tshift = tshift, fig_index = 300, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)
        # os.system('mv *Array_6_pred_wig.png plot.png')

    win_norm = True
    trace_norm = False
    if do_ILAR_pre:
        Zstart_buff = -20 # analysis window start relative to phase arrival
        wind_len    =  18 # analysis window length
        plot_peak = 0.1
        run_compare_pair(repeater = repeater, eq_num1 = eq_num1, eq_num2 = eq_num2, dphase = 'PKIKP', shift_both = shift_both,
                beam_width = beam_width, slow_delta = slow_delta, beam_offset = beam_offset, win_norm = win_norm,trace_norm = trace_norm,
                freq_min = freq_min, freq_max = freq_max, stat_corr = 0, ARRAY = 6, trace_amp = 0.5,
                Zstart_buff = Zstart_buff, wind_len = wind_len, wind_buff = wind_buff, flip = flip,
                tshift = tshift, fig_index = 400, do_interpolate =  True, pair_name = repeater, plot_peak = plot_peak)

    elapsed_time_wc = time.time() - start_time_wc
    print(f'This job took {elapsed_time_wc:.1f} seconds')
    os.system('say "All done"')
