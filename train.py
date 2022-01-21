from browser import document, window, bind, html
 
stA = 'MDGPTRGHGLRKKRRSRSQRDRERRSRGGLGAGAAGGGGAGRTRALSLASSSGSDKEDNGKPPSSAPSRPRPPRRKRRESTSAEEDIIDGFAMTSFVTFEALEKDVALKPQERVEKRQTPLTKKKREALTNGLSFHSKKSRLSHPHHYSSDRENDRNLCQHLGKRKKMPKALRQLKPGQNSCRDSDSESASGESKGFHRSSSRERLSDSSAPSSLGTGYFCDSDSDQEEKASDASSEKLFNTVIVNKDPELGVGTLPEHDSQDAGPIVPKISGLERSQEKSQDCCKEPIFEPVVLKDPCPQVAQPIPQPQTEPQLRAPSPDPDLVQRTEAPPQPPPLSTQPPQGPPEAQLQPAPQPQVQRPPRPQSPTQLLHQNLPPVQAHPSAQSLSQPLSAYNSSSLSLNSLSSSRSSTPAKTQPAPPHISHHPSASPFPLSLPNHSPLHSFTPTLQPPAHSHHPNMFAPPTALPPPPPLTSGSLQVAGHPAGSTYSEQDILRQELNTRFLASQSADRGASLGPPPYLRTEFHQHQHQHQHTHQHTHQHTFTPFPHAIPPTAIMPTPAPPMFDKYPTKVDPFYRHSLFHSYPPAVSGIPPMIPPTGPFGSLQGAFQPKTSNPIDVAARPGTVPHTLLQKDPRLTDPFRPMLRKPGKWCAMHVHIAWQIYHHQQKVKKQMQSDPHKLDFGLKPEFLSRPPGPSLFGAIHHPHDLARPSTLFSAAGAAHPTGTPFGPPPHHSNFLNPAAHLEPFNRPSTFTGLAAVGGNAFGGLGNPSVTPNSMFGHKDGPSVQNFSNPHEPWNRLHRTPPSFPTPPPWLKPGELERSASAAAHDRDRDVDKRDSSVSKDDKERESVEKRHSSHPSPAPVLPVNALGHTRSSTEQIRAHLNTEAREKDKPKERERDHSESRKDLAADEHKAKEGHLPEKDGHGHEGRAAGEEAKQLARVPSPYVRTPVVESARPNSTSSREAEPRKGEPAYENPKKSSEVKVKEERKEDHDLPPEAPQTHRASEPPPPNSSSSVHPGPLASMPMTVGVTGIHPMNSISSLDRTRMMTPFMGISPLPGGERFPYPSFHWDPIRDPLRDPYRELDIHRRDPLGRDFLLRNDPLHRLSTPRLYEADRSFRDREPHDYSHHHHHHHHPLSVDPRREHERGGHLDERERLHMLREDYEHTRLHSVHPASLDGHLPHPSLITPGLPSMHYPRISPTAGNQNGLLNKTPPTAALSAPPPLISTLGGRPVSPRRTTPLSAEIRERPPSHTLKDIEAR'
stA2 = 'PELGVGTLPEHDSQDAGPIVPKISGLERSQEKSQDCCKEPIFEPVVLKDPCPQVAQPIPQPQTEPQLRAPSPDPDLVQRTEAPPQPPPLSTQPPQGPPEAQLQPAPQPQVQRPPRPQSPTQLLHQNLPPVQAHPSAQSLSQPLSAYNSSSLSLNSLSSSRSSTPAKTQPAPPHISHHPSASPFPLSLPNHSPLHSFTPTLQPPAHSHHPNMFAPPTALPPPPPLTSGSLQVAGHPAGSTYSEQDILRQELNTRFLASQSADRGASLGPPPYLRTEFHQHQHQHQHTHQHTHQHTFTPFPH AIPPTAIMPTPAPPMFDKYPTKVDPFYRHSLFHSYPPAVSGIPPMIPPTGPFGSLQGAFQPKLTDPFRPMLRKPGKWCAMHVHIAWQIYHHQQKVKKQMQSDPHKLDFGLKPEFLSRPPG PSLFGAIHHPHDLARPSTLFSAAGAAHPTGTPFGPPPHHSNFLNPAAHLEPFNRPSTFTGLAAVGGNAFGGLGNPSVTPNSMFGHKDGPSVQNFSNPHEPWNRLHRTPPSFPTPPPWLKPGELERSASAAAHDRDRDVDKRDSSVSKDDKERESVEKRHSSHPSPAPVLPVNALGHTRSSTEQIRAHLNTEAREKDKPKERERDHSESRKDLAADEHKAKEGHLPEKDGHGHEGRAAGEEAKQLARVPSPYVRTPVVESARPNSTSSREAEPRKGEPAYENPKKSSEVKVKEERKEDHDLPPEAPQTHRASEPPPPNSSSSVHPGPLASMPMTVGVTGIHPMNSISSLDRTRMMTPFMGISPLPGGERFPYPSFHWDPIRDPLRDPYRELDIHRRDPLGRDFLLRNDPLHRLSTPRLYEADRSFRDREPHDYSHHHHHHHHPLSVDPRREHERGGHLDERERLHMLREDYEHTRLHSVHPASLDGHLPHPSLITPGLPSMHYPRISPTAGNQNGLLNKTPPTAALSAPPPLISTLGGRPVSPRRTTPLSAEIRERPPSHTLKDIEAR'
s = int((len(stA))/40)
s2 = int((len(stA2))/30)
document['seqFstAutismP'].textContent = stA[0:s]
document['seqSndAutismP'].textContent = stA2[0:s2]
#-----------------------------------------------------------------------------------
stADNA = ' CCCCAAAACCCCAAAACCCCTTGCCCTATAAATAGCGAAATAAAGTGGGAAATTATAACTTTATTTTATATTAAAACAAAATCAAATCTTCATCATGAGAGAGGTTATTTCAATTCACGTTGGTCAAGCTGGTATCCAAGTCGGTAACGCATGCTGGGAGCTCTTCTGCCTTGAGCACGGTATTCAACCTGACGGTCAAATGCCATCAAACAAAACCATTGGTGGTGGTGATGATGCCTTCAACACTTTCTTCTCCGAAACTGGAACTGGAAAGCACGTTCCAAGATGCGTCTTCCTCGATTTAGGGCCAACCGTCATTGACGAGGTCAGAACCGGTACCTACAGACAACTCTTCCATCCAGAGCAACTCATCTCAGGAAAGGAGGATGCTGCCAACAACTTCGCCAGAGGTCATTACACCATTGGTAAAGAAATCGTCGATCTCTGCCTCGACGGAATCAGAAAGCTCGCTGATCAATGCACTGGTCTCCAAGGTTTCCTCGTCTTCAACTCAGTTGGTGGTGGTACTGGATCGGGTCTTGGTTCACTCCTCCTCGAGAGACTTTCAGTCGACTATGGTAAGAAATCAAAGCTCGGTTTCACCATCTATCCATCCCCATAAGTCTCAACTGCCGTCGTTGAGCCATACAACTCAGTCCTTTCAACTCACTCACTCCTTGAGCACACTGATGTTGCAGTCATGCTCGACAACGAAGCTGTCTACGATATCTGCAGAAGAAACCTCGACATTGAGAGACCAACCTACACCAACCTTAACAGACTCATCGCTCAAGTTATCTCATCATTGACTGCTTCACTCAGATTCGATGGTGCCCTTAACGTCGATGTTACTGAGTTCCAAACCAACTTGGTCCCATATCCAAGAATCCATTTCATGTTGTCATCCGCCCCAGTCATCTCAGCTGAGAAAGCTTACCACGAGCAACTCTCAGTCGCTGAGATCACCAACTCTGCCTTCGAGCCAGCCTCCATGATGGCCAAGTGCGACCCAAGACACGGTAAATATATGGCTTGCTGCCTTATGTACAGAGGTGATGTCGTCCCCAAGGATGTCAACGCTGCCGTCGCCACCATCA'
stADNA2 = 'TTTTTCCCTTTGTATTTTCCAGAACCACTTGACTGTTTTCCCTACTATACAAGCAACCTATTTTGGAATTATTTCTTTAGCTTTACTCTTTATCAACTATAATCTTTTGGGCAAAGAATTCACCTAAGCTTTTGCTGGTATGAATTATGAGTTCATTCATACTGGCCATTTCAAGTTTAAGCTTTAACTTTCATGCTACCCAGTAAGAAATCTGATAATATGCTATATAATTCATCTTGGTTTACTATATGGTATTCCATTGTACAGCAAATGTATTGTAGGATTCCATTTATATGAAGTACAAAAACAGACAAAAATCATTTATGCTATTAAAAGTCAGGATAGTGGTTACCCTTGATAGGGAGGTGATGACTAAGAGTGGGCACAAATTTGGCAAATAGTATTGTATGAATAAACCACATCTTTTTTTTTTTTTTTTATTTTAAGTTCTGGGATGCATGTCCAGAACGTGCAGGTTTGTTACATAGGTATACATGTGCCATGATGGTTTGCTGCACCTGTCAACCCATCATCTAGGTTTTAAGCTGCGCAAGCATTAGGTATTTGTCCTAATGCTCTCCCTCCCCTTGCCCCCCACCCCTCAACAGGCCCCAGTGTGTGATGTTCCCCTCCCTGTGTCCATGTGTTCTCATTGTTCAACTCCCACTTATGAGTGAGAACATGCGGTGTTTGATTTTCTGTTCCTGTGTTAGTTTGCTGAGAATGATGGCTTCCAGCTTCATTTATGTCCCAAATAAACCACATCTTATTTATCTATTTTACTATTCATAAGCATTGGGTAGTGTCATTAATAATAAAAAAAAACGAAGTATTATATAGGTACTCATATGTACACAAGCAATACTCATATGCATTTATATATACACAGCATTTTAAAGCTTTTTATTGAAGCATAATTTACATACACAGAGATGCATATATCAATAGTACAGCTCAATGAATTTTCACAAACTGAACACACCAGTACCCAGCTCAAGAAACAACACTACCAGCTCCCCAGAAGCCGATTGGTGTTCTCTCCTAGTCATC '
s = int((len(stA))/30)
s2 = int((len(stA2))/20)
document['seqFstAutismD'].textContent = stADNA[20:s]
document['seqSndAutismD'].textContent = stADNA2[20:s2]


storage = window.localStorage
ele = document["Result"] 
eleProtein = document["ResultProtein"]
eleDna = document["ResultDna"]
# ____________________________________________________________________________________________________________________
blosuM62DArray = [
    ['  ', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', ' I ', 'L', 'K','M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'],
    ['A', ' 4 ', '-1 ', '-2 ', '-2 ', ' 0 ', '-1 ', '-1 ', ' 0 ', '-2 ', '-1 ', '-1 ', '-1 ', '-1 ', '-2 ', '-1 ', ' 1 ', ' 0 ',  '-3 ', '-2 ', ' 0 ', '-2 ', '-1 ', ' 0 ', '-4'],
    ['R', '-1 ', ' 5 ', ' 0 ', '-2 ', '-3 ', ' 1 ', ' 0 ', '-2 ', ' 0 ', '-3 ', '-2 ', ' 2 ',    '-1 ', '-3 ', '-2 ', '-1 ', '-1 ', '-3 ', '-2 ', '-3 ', '-1 ', ' 0 ', '-1 ', '-4'],
    ['N', '-2 ', ' 0 ', ' 6 ', ' 1 ', '-3 ', ' 0 ', ' 0 ', ' 0 ', ' 1 ', '-3 ', '-3 ', ' 0 ',     '-2 ', '-3 ', '-2 ', ' 1 ', ' 0 ', '-4 ', '-2 ', '-3 ', ' 3 ', ' 0 ', '-1 ', '-4'],
    ['D', '-2 ', '-2 ', ' 1 ', ' 6 ', '-3 ', ' 0 ', ' 2 ', '-1 ', '-1 ', '-3 ', '-4 ', '-1 ',     '-3 ', '-3 ', '-1 ', ' 0 ', '-1 ', '-4 ', '-3 ', '-3 ', ' 4 ', ' 1 ', '-1 ', '-4'],
    ['C', ' 0 ', '-3 ', '-3 ', '-3 ', ' 9 ', '-3 ', '-4 ', '-3 ', '-3 ', '-1 ', '-1 ', '-3 ',     '-1 ', '-2 ', '-3 ', '-1 ', '-1 ', '-2 ', '-2 ', '-1 ', '-3 ', '-3 ', '-2 ', '-4'],
    ['Q', '-1 ', ' 1 ', ' 0 ', ' 0 ', '-3 ', ' 5 ', ' 2 ', '-2 ', ' 0 ', '-3 ', '-2 ', ' 1 ',    ' 0 ', '-3 ', '-1 ', ' 0 ', '-1 ', '-2 ', '-1 ', '-2 ', ' 0 ', ' 3 ', '-1 ', '-4'],
    ['E', '-1 ', ' 0 ', ' 0 ', ' 2 ', '-4 ', ' 2 ', ' 5 ', '-2 ', ' 0 ', '-3 ', '-3 ', ' 1 ',    '-2 ', '-3 ', '-1 ', ' 0 ', '-1 ', '-3 ', '-2 ', '-2 ', ' 1 ', ' 4 ', '-1 ', '-4'],
    ['G', ' 0 ', '-2 ', ' 0 ', '-1 ', '-3 ', '-2 ', '-2 ', ' 6 ', '-2 ', '-4 ', '-4 ', '-2 ',    '-3 ', '-3 ', '-2 ', ' 0 ', '-2 ', '-2 ', '-3 ', '-3 ', '-1 ', '-2 ', '-1 ', '-4'],
    ['H', '-2 ', ' 0 ', ' 1 ', '-1 ', '-3 ', ' 0 ', ' 0 ', '-2 ', ' 8 ', '-3 ', '-3 ', '-1 ',      '-2 ', '-1 ', '-2 ', '-1 ', '-2 ', '-2 ', ' 2 ', '-3 ', ' 0 ', ' 0 ', '-1 ', '-4'],
    ['I', '-1 ', '-3 ', '-3 ', '-3 ', '-1 ', '-3 ', '-3 ', '-4 ', '-3 ', ' 4 ', ' 2 ', '-3 ',    ' 1 ', ' 0 ', '-3 ', '-2 ', '-1 ',  '-3 ', '-1 ', ' 3 ', '-3 ', '-3 ', '-1 ', '-4'],
    ['L', '-1 ', '-2 ', '-3 ', '-4 ', '-1 ', '-2 ', '-3 ', '-4 ', '-3 ', ' 2 ', ' 4 ', '-2 ',       ' 2 ', ' 0 ', '-3 ', '-2 ', '-1 ', '-2 ', '-1 ', ' 1 ', '-4 ', '-3 ', '-1 ', '-4'],
    ['K', '-1 ', ' 2 ', ' 0 ', '-1 ', '-3 ', ' 1 ', ' 1 ', '-2 ', '-1 ', '-3 ', '-2 ', ' 5 ',     '-1 ', '-3 ', '-1 ', ' 0 ', '-1 ', '-3 ', '-2 ', '-2 ', ' 0 ', ' 1 ', '-1 ', '-4'],
    ['M', '-1 ', '-1 ', '-2 ', '-3 ', '-1 ', ' 0 ', '-2 ', '-3 ', '-2 ', ' 1 ', ' 2 ', '-1 ',     ' 5 ', ' 0 ', '-2 ', '-1 ', '-1 ',      '-1 ', '-1 ', ' 1 ', '-3 ', '-1 ', '-1 ', '-4'],
    ['F', '-2 ', '-3 ', '-3 ', '-3 ', '-2 ', '-3 ', '-3 ', '-3 ', '-1 ', ' 0 ', ' 0 ', '-3 ',     ' 0 ', ' 6 ', '-4 ', '-2 ', '-2 ',      ' 1 ', ' 3 ', '-1 ', '-3 ', '-3 ', '-1 ', '-4'],
    ['P', '-1 ', '-2 ', '-2 ', '-1 ', '-3 ', '-1 ', '-1 ', '-2 ', '-2 ', '-3 ', '-3 ', '-1 ',
     '-2 ', '-4 ', ' 7 ', '-1 ', '-1 ',      '-4 ', '-3 ', '-2 ', '-2 ', '-1 ', '-2 ', '-4'],
    ['S', ' 1 ', '-1 ', ' 1 ', ' 0 ', '-1 ', ' 0 ', ' 0 ', ' 0 ', '-1 ', '-2 ', '-2 ', ' 0 ',
     '-1 ', '-2 ', '-1 ', ' 4 ', ' 1 ',     '-3 ', '-2 ', '-2 ', ' 0 ', ' 0 ', ' 0 ', '-4'],
    ['T', ' 0 ', '-1 ', ' 0 ', '-1 ', '-1 ', '-1 ', '-1 ', '-2 ', '-2 ', '-1 ', '-1 ', '-1 ',
     '-1 ', '-2 ', '-1 ', ' 1 ', ' 5 ',      '-2 ', '-2 ', ' 0 ', '-1 ', '-1 ', ' 0 ', '-4'],
    ['W ', '-3 ', '-3 ', '-4 ', '-4 ', '-2 ', '-2 ', '-3 ', '-2 ', '-2 ', '-3 ', '-2 ', '-3 ',
           '-1 ', ' 1 ', '-4 ', '-3 ', '-2 ',       '11 ', ' 2 ', '-3 ', '-4 ', '-3 ', '-2 ', '-4'],
    ['Y ', '-2 ', '-2 ', '-2 ', '-3 ', '-2 ', '-1 ', '-2 ', '-3 ', ' 2 ', '-1 ', '-1 ', '-2 ',
           '-1 ', ' 3 ', '-3 ', '-2 ', '-2 ',     ' 2 ', ' 7 ', '-1 ', '-3 ', '-2 ', '-1 ', '-4'],
    ['V ', ' 0 ', '-3 ', '-3 ', '-3 ', '-1 ', '-2 ', '-2 ', '-3 ', '-3 ', ' 3 ', ' 1 ', '-2 ',
           ' 1 ', '-1 ', '-2 ', '-2 ', ' 0 ',    '-3 ', '-1 ', ' 4 ', '-3 ', '-2 ', '-1 ', '-4'],
    ['B ', '-2 ', '-1 ', ' 3 ', ' 4 ', '-3 ', ' 0 ', ' 1 ', '-1 ', ' 0 ', '-3 ', '-4 ', ' 0 ',
           '-3 ', '-3 ', '-2 ', ' 0 ', '-1 ',     '-4 ', '-3 ', '-3 ', ' 4 ', ' 1 ', '-1 ', '-4'],
    ['Z ', '-1 ', ' 0 ', ' 0 ', ' 1 ', '-3 ', ' 3 ', ' 4 ', '-2 ', ' 0 ', '-3 ', '-3 ', ' 1 ',
           '-1 ', '-3 ', '-1 ', ' 0 ', '-1 ',        '-3 ', '-2 ', '-2 ', ' 1 ', ' 4 ', '-1 ', '-4'],
    ['X ', ' 0 ', '-1 ', '-1 ', '-1 ', '-2 ', '-1 ', '-1 ', '-1 ',  '-1 ', '-1 ', '-1 ', '-1 ',
           '-1 ', '-1 ', '-2 ', ' 0 ', ' 0 ',    '-2 ', '-1 ', '-1 ', '-1 ', '-1 ', '-1 ', '-4'],
    ['* ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ',  '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ',     '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', '-4 ', ' 1']]
# ____________________________________________________________________________________________________________________
def globalAlignment(e):
    global ele
    #ele.style.visibility = "visible"
    ele.style.display = "block"
   # _____________________________
    
    print(s,s2)
    sequenceFst = document['seq1-input'] .value.upper()
    sequenceSnd = document['seq2-input'].value.upper()
    
    # _____________________________

    match_bonus = int(document['match-input'].value)
    mismatch_penalty = int(document['misMatch-input'].value)
    gap_penalty = int(document['gap-input'].value)

    # _____________________________
    document['displaySeq1'].textContent = sequenceFst
    document['displaySeq2'].textContent = sequenceSnd
    storage.setItem('displaySeq1', sequenceFst)
    storage.setItem('displaySeq2', sequenceSnd)

    #  storage.setItem('seq2', seq2)
    n_rows = len("-"+sequenceFst)
    print(len(sequenceFst))
    print(n_rows)
    n_columns = len("-"+sequenceSnd)
   
    scoring_array = [([0]*n_columns) for i in range(n_rows)]
    traceback_array = [(['-']*n_columns) for i in range(n_rows)]

    # _________________________________________________________
    # Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2193"
    #right_arrow = "\u2190"
    #down_arrow = "\u2191"
    left_arrow = "\u2192"
    #down_right_arrow = "\u2196"
    up_left_arrow = "\u2198"
    arrow = "-"
    # _________________________________________________________
    # _________________________________________________________

    for row in range(n_rows):
        for col in range(n_columns):
            if row == 0 and col == 0:
                # We're in the upper right corner
                score = 0
                arrow = "-"
            elif row == 0:             
                previous_score = scoring_array[row][col - 1]            
                score = previous_score + gap_penalty
                arrow = left_arrow
            elif col == 0:            
                previous_score = scoring_array[row - 1][col]
                score = previous_score + gap_penalty
                arrow = up_arrow
            else:            
                cell_to_the_left = scoring_array[row][col-1]
                from_left_score = cell_to_the_left + gap_penalty

                above_cell = scoring_array[row-1][col]
                from_above_score = above_cell + gap_penalty
                
                diagonal_left_cell = scoring_array[row-1][col-1]

                if sequenceFst[row-1] == sequenceSnd[col-1]:
                    diagonal_left_cell_score = diagonal_left_cell + match_bonus
                else:
                    diagonal_left_cell_score = diagonal_left_cell + mismatch_penalty

                score = max([from_left_score, from_above_score,
                            diagonal_left_cell_score])
                # take the max
                scoring_array[row][col] = score

                # using Unicode arrows
                if score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow           
            traceback_array[row][col] = arrow
            
            document['tracebackArray'].textContent = [ ('\n', traceback_array[i]) for i in range(len(traceback_array))]
            document['scoreArray'].textContent = [('\n', scoring_array[i]) for i in range(len(scoring_array))]

            print(score)

            scoring_array[row][col] = score
            print(scoring_array)
    traceback_alignment(traceback_array, sequenceFst, sequenceSnd)
    document['scoredAlignment'].textContent = '\n', traceback_alignment( traceback_array, sequenceFst, sequenceSnd)[0], '\n',traceback_alignment( traceback_array, sequenceFst, sequenceSnd)[1], '\n',traceback_alignment( traceback_array, sequenceFst, sequenceSnd)[2]
# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________


def globalAlignmentProtein(e):
    global ele
    ele.style.display = "block"

    # #_____________________________
    
    sequenceFst = document['seq1-inputPP'] .value.upper()
    sequenceSnd = document['seq2-inputPP'].value.upper()
    # _____________________________
    gap_penalty = int(document['gap-inputPP'].value)

    # #_____________________________
    document['displaySeq1'].textContent = sequenceFst
    document['displaySeq2'].textContent = sequenceSnd
    storage.setItem('displaySeq1', sequenceFst)
    storage.setItem('displaySeq2', sequenceSnd)

    #  storage.setItem('seq2', seq2)
    n_rows = len("-"+sequenceFst)
   # print(len(sequenceFst))
   # print(  n_rows)
    n_columns = len("-"+sequenceSnd)

    scoring_array = [([0]*n_columns) for i in range(n_rows)]
    traceback_array = [(['-']*n_columns) for i in range(n_rows)]

    # _________________________________________________________
    # Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2193"
    #right_arrow = "\u2190"
    #down_arrow = "\u2191"
    left_arrow = "\u2192"
    #down_right_arrow = "\u2196"
    up_left_arrow = "\u2198"
    arrow = "-"
    # _________________________________________________________
    # _________________________________________________________

    for row in range(n_rows):
        for col in range(n_columns):
            if row == 0 and col == 0:
                # We're in the upper right corner
                score = 0
                arrow = "-"
            elif row == 0:
                
                previous_score = scoring_array[row][col - 1]
 
                score = previous_score + gap_penalty
                arrow = left_arrow
            elif col == 0:
                
                previous_score = scoring_array[row - 1][col]
                score = previous_score + gap_penalty
                arrow = up_arrow
            else:
                cell_to_the_left = scoring_array[row][col-1]
                from_left_score = cell_to_the_left + gap_penalty

                above_cell = scoring_array[row-1][col]
                from_above_score = above_cell + gap_penalty
               #_______________________________________________                              
                diagonal_left_cell = scoring_array[row-1][col-1]
                
                
                print("Hy global row",row)

                print(sequenceFst[row-1], '====equenceFst Align====', row, ' \n',
                        sequenceSnd[col-1], '===sequenceSnd Align====', col, ' \n')
               
                check = sequenceFst[row-1] in blosuM62DArray[0][0:]
                check2 = sequenceSnd[col-1] in blosuM62DArray[0:][0]
                if (check == True & check2 == True):
                    match_bonus = int(
                        blosuM62DArray
                        [blosuM62DArray[0][0:].index(sequenceFst[row-1])]
                        [blosuM62DArray[0:][0].index(sequenceSnd[col-1])])
                print("match_bonus",match_bonus)
               

                #_______________________________________________
                diagonal_left_cell_score = diagonal_left_cell+  match_bonus

                score = max([from_left_score, from_above_score,diagonal_left_cell_score])
                print("score",score)
                print("from_left_score",from_left_score)
                print("from_above_score",from_above_score)
                print("diagonal_left_cell_score",diagonal_left_cell_score)
                print("scoring_array[rowA][colA]",scoring_array[row-1][col-1])

                # take the max
                scoring_array[row][col] = score
                print("scoring_array[row][col]",scoring_array[row][col])

                if score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow

            traceback_array[row][col] = arrow
            # print(traceback_array)
            document['tracebackArray'].textContent = [
                ('\n', traceback_array[i]) for i in range(len(traceback_array))]
            document['scoreArray'].textContent = [('\n', scoring_array[i] ) for i in range(len(scoring_array))]

           # print(score)

            scoring_array[row][col] = score
            # print(scoring_array)
    traceback_alignment(traceback_array, sequenceFst, sequenceSnd)

    document['scoredAlignment'].textContent = '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[0], '\n', traceback_alignment(
        traceback_array, sequenceFst, sequenceSnd)[1], '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[2]

# ____________________________________________________________________________________________________________________


def traceback_alignment(traceback_array, seq1, seq2, up_arrow="\u2193", left_arrow="\u2192", up_left_arrow="\u2198", stop="-"):
    n_rows = len(seq1) + 1  # need an extra row up top
    n_columns = len(seq2) + 1  # need an extra row up top

    row = len(seq1)
    col = len(seq2)
    arrow = traceback_array[row][col]
    aligned_seq1 = ""
    aligned_seq2 = ""
    alignment_indicator = ""
    while arrow is not "-":
        print("Currently on row:", row)
        print("Currently on col:", col)
        arrow = traceback_array[row][col]
        print("Arrow:", arrow)

        if arrow == up_arrow:
            print("insert indel into top sequence")      
            aligned_seq2 = "-"+aligned_seq2
            aligned_seq1 = seq1[row-1] + aligned_seq1
            alignment_indicator = " "+alignment_indicator
            row -= 1

        elif arrow == up_left_arrow:
            print("match or mismatch")         
            seq1_character = seq1[row-1]
            seq2_character = seq2[col-1]
            aligned_seq1 = seq1[row-1] + aligned_seq1
            aligned_seq2 = seq2[col-1] + aligned_seq2
            if seq1_character == seq2_character:
                alignment_indicator = "|"+alignment_indicator
            else:
                alignment_indicator = " "+alignment_indicator
            row -= 1
            col -= 1

        elif arrow == left_arrow:
            print("Insert indel into left sequence")
            aligned_seq1 = "-"+aligned_seq1
            aligned_seq2 = seq2[col-1] + aligned_seq2
            alignment_indicator = " "+alignment_indicator
            col -= 1

        elif arrow == stop:
            break
        else:
            raise ValueError(
                f"Traceback array entry at {row},{col}: {arrow} is not recognized as an up arrow ({up_arrow}),left_arrow ({left_arrow}), up_left_arrow ({up_left_arrow}), or a stop ({stop}).")
        print(aligned_seq1.replace(',',''))
        print(alignment_indicator.replace(',', ''))
        print(aligned_seq2.replace(',', ''))

    return aligned_seq1, alignment_indicator,aligned_seq2
# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________


def globalAutismAlignmentDNA(e):
    global eleDna
    eleDna.style.display = "block"
   # _____________________________

    sequenceFst = document['seqFstAutismD'].innerText
    sequenceSnd = document['seqSndAutismD'].innerText
    print('testing..', sequenceFst, sequenceSnd)
    
    
    
    match_bonus = int(document['match-inputD'].value)
    mismatch_penalty = int(document['misMatch-inputD'].value)
    gap_penalty = int(document['gap-inputD'].value)

    # _____________________________
    document['displaySeq1D'].textContent = sequenceFst
    document['displaySeq2D'].textContent = sequenceSnd
    storage.setItem('displaySeq1D', sequenceFst)
    storage.setItem('displaySeq2D', sequenceSnd)

    #  storage.setItem('seq2', seq2)
    n_rows = len("-"+sequenceFst)
    print(len(sequenceFst))
    print(n_rows)
    n_columns = len("-"+sequenceSnd)

    scoring_array = [([0]*n_columns) for i in range(n_rows)]
    traceback_array = [(['-']*n_columns) for i in range(n_rows)]

    # _________________________________________________________
    # Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2193"
    #right_arrow = "\u2190"
    #down_arrow = "\u2191"
    left_arrow = "\u2192"
    #down_right_arrow = "\u2196"
    up_left_arrow = "\u2198"
    arrow = "-"
    # _________________________________________________________
    # _________________________________________________________

    for row in range(n_rows):
        for col in range(n_columns):
            if row == 0 and col == 0:
                # We're in the upper right corner
                score = 0
                arrow = "-"
            elif row == 0:
                previous_score = scoring_array[row][col - 1]
                score = previous_score + gap_penalty
                arrow = left_arrow
            elif col == 0:
                previous_score = scoring_array[row - 1][col]
                score = previous_score + gap_penalty
                arrow = up_arrow
            else:
                cell_to_the_left = scoring_array[row][col-1]
                from_left_score = cell_to_the_left + gap_penalty

                above_cell = scoring_array[row-1][col]
                from_above_score = above_cell + gap_penalty

                diagonal_left_cell = scoring_array[row-1][col-1]

                if sequenceFst[row-1] == sequenceSnd[col-1]:
                    diagonal_left_cell_score = diagonal_left_cell + match_bonus
                else:
                    diagonal_left_cell_score = diagonal_left_cell + mismatch_penalty

                score = max([from_left_score, from_above_score,
                            diagonal_left_cell_score])
                # take the max
                scoring_array[row][col] = score

                # using Unicode arrows
                if score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow

            traceback_array[row][col] = arrow

            document['tracebackArrayD'].textContent = [
                ('\n', traceback_array[i]) for i in range(len(traceback_array))]
            document['scoreArrayD'].textContent = [
                ('\n', scoring_array[i]) for i in range(len(scoring_array))]

            print(score)

            scoring_array[row][col] = score
            print(scoring_array)
    traceback_alignment(traceback_array, sequenceFst, sequenceSnd)
    document['scoredAlignmentD'].textContent = '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[0], '\n', traceback_alignment(
        traceback_array, sequenceFst, sequenceSnd)[1], '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[2]


def globalAutismAlignmentProtein(e):
    global eleProtein
    eleProtein.style.display = "block"

    # #_____________________________

    sequenceFst = document['seqFstAutismP'].innerText
    sequenceSnd = document['seqSndAutismP'].innerText
    print('oooooooooooooooo',sequenceFst, sequenceSnd)

    document['displaySeq1P'].textContent = sequenceFst
    document['displaySeq2P'].textContent = sequenceSnd
    storage.setItem('displaySeq1P', sequenceFst)
    storage.setItem('displaySeq2P', sequenceSnd)

    # _____________________________

    gap_penalty = int(document['gap-inputP'].value)

    # #_____________________________
    #  storage.setItem('seq2', seq2)
    n_rows = len("-"+sequenceFst)
   # print(len(sequenceFst))
   # print(  n_rows)
    n_columns = len("-"+sequenceSnd)

    scoring_array = [([0]*n_columns) for i in range(n_rows)]
    traceback_array = [(['-']*n_columns) for i in range(n_rows)]

    # _________________________________________________________
    # Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2193"
    #right_arrow = "\u2190"
    #down_arrow = "\u2191"
    left_arrow = "\u2192"
    #down_right_arrow = "\u2196"
    up_left_arrow = "\u2198"
    arrow = "-"
    # _________________________________________________________
    # _________________________________________________________

    for row in range(n_rows):
        for col in range(n_columns):
            if row == 0 and col == 0:
                # We're in the upper right corner
                score = 0
                arrow = "-"
            elif row == 0:
                # We're on the first row
                # but NOT in the corner

                # Look up the score of the previous cell (to the left) in the score array\
                previous_score = scoring_array[row][col - 1]
                # add the gap penalty to it's score
                score = previous_score + gap_penalty
                arrow = left_arrow
            elif col == 0:
                # We're on the first column but not in the first row
                previous_score = scoring_array[row - 1][col]
                score = previous_score + gap_penalty
                arrow = up_arrow
            else:
                cell_to_the_left = scoring_array[row][col-1]
                from_left_score = cell_to_the_left + gap_penalty

                above_cell = scoring_array[row-1][col]
                from_above_score = above_cell + gap_penalty
               #_______________________________________________

                diagonal_left_cell = scoring_array[row-1][col-1]
                print("Hy global row", row)

                print(sequenceFst[row-1], '====equenceFst Align====', row, ' \n',
                      sequenceSnd[col-1], '===sequenceSnd Align====', col, ' \n')

                check = sequenceFst[row-1] in blosuM62DArray[0][0:]
                check2 = sequenceSnd[col-1] in blosuM62DArray[0:][0]
                if (check == True & check2 == True):
                    match_bonus = int(
                        blosuM62DArray
                        [blosuM62DArray[0][0:].index(sequenceFst[row-1])]
                        [blosuM62DArray[0:][0].index(sequenceSnd[col-1])])
                print("match_bonus", match_bonus)

                #_______________________________________________
                diagonal_left_cell_score = diagonal_left_cell + match_bonus

                score = max([from_left_score, from_above_score,
                            diagonal_left_cell_score])
                print("score", score)
                print("from_left_score", from_left_score)
                print("from_above_score", from_above_score)
                print("diagonal_left_cell_score", diagonal_left_cell_score)
                print("scoring_array[rowA][colA]", scoring_array[row-1][col-1])

                # take the max
                scoring_array[row][col] = score
                print("scoring_array[row][col]", scoring_array[row][col])

                if score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow

            traceback_array[row][col] = arrow
            # print(traceback_array)
            document['tracebackArrayP'].textContent = [
                ('\n', traceback_array[i]) for i in range(len(traceback_array))]
            document['scoreArrayP'].textContent = [
                ('\n', scoring_array[i]) for i in range(len(scoring_array))]

           # print(score)

            scoring_array[row][col] = score
            # print(scoring_array)
    traceback_alignment(traceback_array, sequenceFst, sequenceSnd)

    document['scoredAlignmentP'].textContent = '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[0], '\n', traceback_alignment(
        traceback_array, sequenceFst, sequenceSnd)[1], '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[2]

# ____________________________________________________________________________________________________________________


def localAlignment(e):
    global ele
    ele.style.visibility = "visible"

    # _____________________________

    sequenceFst = document['seq1-input'] .value
    sequenceSnd = document['seq2-input'].value

    # _____________________________

    match_bonus = int(document['match-input'].value)
    mismatch_penalty = int(document['misMatch-input'].value)
    gap_penalty = int(document['gap-input'].value)

    # _____________________________
    document['displaySeq1'].textContent = sequenceFst
    document['displaySeq2'].textContent = sequenceSnd
    storage.setItem('displaySeq1', sequenceFst)
    storage.setItem('displaySeq2', sequenceSnd)

    #  storage.setItem('seq2', seq2)
    n_rows = len("-"+sequenceFst)
    print(len(sequenceFst))
    print(n_rows)
    n_columns = len("-"+sequenceSnd)
    row_labels = [label for label in "-"+sequenceFst]
    column_labels = [label for label in "-"+sequenceSnd]
    #   SetMatrix null
    scoring_array = [([0]*n_columns) for i in range(n_rows)]
    traceback_array = [(['-']*n_columns) for i in range(n_rows)]

    # _________________________________________________________
    # Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2191"
    #right_arrow = "\u2192"
    #down_arrow = "\u2193"
    left_arrow = "\u2190"
    #down_right_arrow = "\u2198"
    up_left_arrow = "\u2196"
    arrow = "-"
    # _________________________________________________________
    # _________________________________________________________

    for row in range(n_rows):
        for col in range(n_columns):
            if row == 0 and col == 0:
                # We're in the upper right corner
                score = 0
                arrow = "-"
            elif row == 0:
                score = 0
                arrow = '-'

            elif col == 0:
                score = 0
                arrow = '-'
            else:

                cell_to_the_left = scoring_array[row][col-1]
                from_left_score = cell_to_the_left + gap_penalty
                above_cell = scoring_array[row-1][col]
                from_above_score = above_cell + gap_penalty

                diagonal_left_cell = scoring_array[row-1][col-1]

                if sequenceFst[row-1] == sequenceSnd[col-1]:
                    diagonal_left_cell_score = diagonal_left_cell + match_bonus
                else:
                    diagonal_left_cell_score = diagonal_left_cell + mismatch_penalty

                score = max([from_left_score, from_above_score,
                            diagonal_left_cell_score])
                # take the max
                if score < 0:
                    score = 0
                scoring_array[row][col] = score
                # print(score)
                # make note of which cell was the max in the traceback array
                # using Unicode arrows
                if score == 0:
                    arrow = '-'
                elif score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow
            # print(arrow)
            # print(traceback_array)
            # print(traceback_array[1][1])

            traceback_array[row][col] = arrow
            # print(traceback_array)
            document['tracebackArray'].textContent = [
                ('\n', traceback_array[i]) for i in range(len(traceback_array))]
            document['scoreArray'].textContent = [
                ('\n', scoring_array[i]) for i in range(len(scoring_array))]

            print(score)

            scoring_array[row][col] = score
            print(scoring_array)
    maxScore = max(max(x) for x in scoring_array)
 
    for i in range(len(scoring_array)):
        for j in range(len(scoring_array)):
            if(scoring_array[i][j] == maxScore):
                print(scoring_array[i][j], i, j)
                rowMaxAlign = i
                colMaxAlign = j
                strAlignment = traceback_local_alignment(
                    traceback_array, sequenceFst, sequenceSnd, rowMaxAlign, colMaxAlign)
                break
 
    document['scoredAlignment'].textContent = '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[0], '\n', traceback_alignment(
        traceback_array, sequenceFst, sequenceSnd)[1], '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[2]


# ____________________________________________________________________________________________________________________


def traceback_local_alignment(traceback_array, seq1, seq2, rowMaxAlign, colMaxAlign, up_arrow="\u2191", left_arrow="\u2190", up_left_arrow="\u2196", stop="-"):
    n_rows = len(seq1) + 1  # need an extra row up top
    n_columns = len(seq2) + 1  # need an extra row up top

    print("Hy in traceback_local_alignment :")
    row = rowMaxAlign
    col = colMaxAlign
    arrow = traceback_array[row][col]
    print(traceback_array[row][col], "max", row, col)
    aligned_seq1 = ""
    aligned_seq2 = ""
    alignment_indicator = ""
    while arrow is not "-":
        print("Currently on row:", row)
        print("Currently on col:", col)
        arrow = traceback_array[row][col]
        print("Arrow:", arrow)

        if arrow == up_arrow:
            print("insert indel into top sequence")

            aligned_seq2 = "-"+aligned_seq2
            aligned_seq1 = seq1[row-1] + aligned_seq1
            alignment_indicator = " "+alignment_indicator
            row -= 1

        elif arrow == up_left_arrow:
            print("match or mismatch")
            # Note that we look up the row-1 and col-1 indexes
            # because there is an extra "-" character at the
            # start of each sequence
            seq1_character = seq1[row-1]
            seq2_character = seq2[col-1]
            aligned_seq1 = seq1[row-1] + aligned_seq1
            aligned_seq2 = seq2[col-1] + aligned_seq2
            if seq1_character == seq2_character:
                alignment_indicator = "|"+alignment_indicator
            else:
                alignment_indicator = " "+alignment_indicator
            row -= 1
            col -= 1

        elif arrow == left_arrow:
            print("Insert indel into left sequence")
            aligned_seq1 = "-"+aligned_seq1
            aligned_seq2 = seq2[col-1] + aligned_seq2
            alignment_indicator = " "+alignment_indicator
            col -= 1

        elif arrow == stop:
            break
        else:
            raise ValueError(
                f"Traceback array entry at {row},{col}: {arrow} is not recognized as an up arrow ({up_arrow}),left_arrow ({left_arrow}), up_left_arrow ({up_left_arrow}), or a stop ({stop}).")
        # print(traceback_array,-row,-col,traceback_array[-row,-col])
        print(aligned_seq1)
        print(alignment_indicator)
        print(aligned_seq2)

    return aligned_seq1, aligned_seq2, alignment_indicator
# ____________________________________________________________________________________________________________________


def localAlignment(e):
    global ele
    ele.style.display = "block"

    # _____________________________

    sequenceFst = document['seq1-input'] .value
    sequenceSnd = document['seq2-input'].value

    # _____________________________

    match_bonus = int(document['match-input'].value)
    mismatch_penalty = int(document['misMatch-input'].value)
    gap_penalty = int(document['gap-input'].value)

    # _____________________________
    document['displaySeq1'].textContent = sequenceFst
    document['displaySeq2'].textContent = sequenceSnd
    storage.setItem('displaySeq1', sequenceFst)
    storage.setItem('displaySeq2', sequenceSnd)
    print(sequenceFst, sequenceSnd)
    print(len(sequenceFst), len(sequenceSnd))

    if(len(sequenceFst) > len(sequenceSnd)):

        sequenceSnd += '-'*(len(sequenceFst) - len(sequenceSnd))

    elif(len(sequenceSnd) > len(sequenceFst)):
        sequenceFst = sequenceFst + '-'*(len(sequenceSnd) - len(sequenceFst))

    print(sequenceFst, sequenceSnd)
    print(len(sequenceFst), len(sequenceSnd))

    #  storage.setItem('seq2', seq2)
    n_rows = len("-"+sequenceFst)
    print(len(sequenceFst))
    print(n_rows)
    n_columns = len("-"+sequenceSnd)
    row_labels = [label for label in "-"+sequenceFst]
    column_labels = [label for label in "-"+sequenceSnd]
    #   SetMatrix null
    scoring_array = [([0]*n_columns) for i in range(n_rows)]
    traceback_array = [(['-']*n_columns) for i in range(n_rows)]

    # _________________________________________________________
    # Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2191"
    #right_arrow = "\u2192"
    #down_arrow = "\u2193"
    left_arrow = "\u2190"
    #down_right_arrow = "\u2198"
    up_left_arrow = "\u2196"
    arrow = "-"
    # _________________________________________________________
    # _________________________________________________________

    for row in range(n_rows):
        for col in range(n_columns):
            if row == 0 and col == 0:
                # We're in the upper right corner
                score = 0
                arrow = "-"
            elif row == 0:
                score = 0
                arrow = '-'

            elif col == 0:
                score = 0
                arrow = '-'
            else:

                cell_to_the_left = scoring_array[row][col-1]
                from_left_score = cell_to_the_left + gap_penalty
                above_cell = scoring_array[row-1][col]
                from_above_score = above_cell + gap_penalty

                diagonal_left_cell = scoring_array[row-1][col-1]

                if sequenceFst[row-1] == sequenceSnd[col-1]:
                    diagonal_left_cell_score = diagonal_left_cell + match_bonus
                else:
                    diagonal_left_cell_score = diagonal_left_cell + mismatch_penalty

                score = max([from_left_score, from_above_score,
                            diagonal_left_cell_score])
                # take the max
                if score < 0:
                    score = 0
                scoring_array[row][col] = score
                # print(score)
                # make note of which cell was the max in the traceback array
                # using Unicode arrows
                if score == 0:
                    arrow = '-'
                elif score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow
            # print(arrow)
            # print(traceback_array)
            # print(traceback_array[1][1])

            traceback_array[row][col] = arrow
            # print(traceback_array)
            document['tracebackArray'].textContent = [
                ('\n', traceback_array[i]) for i in range(len(traceback_array))]
            document['scoreArray'].textContent = [
                ('\n', scoring_array[i]) for i in range(len(scoring_array))]

            print(score)

            scoring_array[row][col] = score
            print(scoring_array)
    maxScore = max(max(x) for x in scoring_array)
   

    for i in range(len(scoring_array)):
        for j in range(len(scoring_array)):
            if(scoring_array[i][j] == maxScore):
                print(scoring_array[i][j], i, j)
                rowMaxAlign = i
                colMaxAlign = j
                strAlignment = traceback_local_alignment(
                    traceback_array, sequenceFst, sequenceSnd, rowMaxAlign, colMaxAlign)
                break
 
    document['scoredAlignment'].textContent = '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[0], '\n', traceback_alignment(
        traceback_array, sequenceFst, sequenceSnd)[1], '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[2]

# ____________________________________________________________________________________________________________________


def localProteinAlignment(e):
    global ele
    ele.style.display = "block"

    # _____________________________

    sequenceFst = document['seq1-inputPP'] .value.upper()
    sequenceSnd = document['seq2-inputPP'].value.upper()

    # _____________________________

    gap_penalty = int(document['gap-inputPP'].value)

    # _____________________________

    # _____________________________
    document['displaySeq1'].textContent = sequenceFst
    document['displaySeq2'].textContent = sequenceSnd
    storage.setItem('displaySeq1', sequenceFst)
    storage.setItem('displaySeq2', sequenceSnd)
    print(sequenceFst, sequenceSnd)
    print(len(sequenceFst), len(sequenceSnd))

    if(len(sequenceFst) > len(sequenceSnd)):

        sequenceSnd += '-'*(len(sequenceFst) - len(sequenceSnd))

    elif(len(sequenceSnd) > len(sequenceFst)):
        sequenceFst = sequenceFst + '-'*(len(sequenceSnd) - len(sequenceFst))

    print(sequenceFst, sequenceSnd)
    print(len(sequenceFst), len(sequenceSnd))

    #  storage.setItem('seq2', seq2)
    n_rows = len("-"+sequenceFst)
    print(len(sequenceFst))
    print(n_rows)
    n_columns = len("-"+sequenceSnd)
    row_labels = [label for label in "-"+sequenceFst]
    column_labels = [label for label in "-"+sequenceSnd]
    #   SetMatrix null
    scoring_array = [([0]*n_columns) for i in range(n_rows)]
    traceback_array = [(['-']*n_columns) for i in range(n_rows)]

    # _________________________________________________________
    # Define Unicode arrows we'll use in the traceback array
    up_arrow = "\u2191"
    #right_arrow = "\u2192"
    #down_arrow = "\u2193"
    left_arrow = "\u2190"
    #down_right_arrow = "\u2198"
    up_left_arrow = "\u2196"
    arrow = "-"
    # _________________________________________________________
    # _________________________________________________________

    for row in range(n_rows):
        for col in range(n_columns):
            if row == 0 and col == 0:
                # We're in the upper right corner
                score = 0
                arrow = "-"
            elif row == 0:
                score = 0
                arrow = '-'

            elif col == 0:
                score = 0
                arrow = '-'
            else:

                cell_to_the_left = scoring_array[row][col-1]
                from_left_score = cell_to_the_left + gap_penalty
                above_cell = scoring_array[row-1][col]
                from_above_score = above_cell + gap_penalty

                diagonal_left_cell = scoring_array[row-1][col-1]

                check = sequenceFst[row-1] in blosuM62DArray[0][0:]
                check2 = sequenceSnd[col-1] in blosuM62DArray[0:][0]
                if (check == True & check2 == True):
                    match_bonus = int(
                        blosuM62DArray
                        [blosuM62DArray[0][0:].index(sequenceFst[row-1])]
                        [blosuM62DArray[0:][0].index(sequenceSnd[col-1])])
                    print("match_bonus", match_bonus)
                    #  if sequenceFst[row-1] == sequenceSnd[col-1]:
                    #print(diagonal_left_cell, '*******zffffftt*')

                    #_______________________________________________
                diagonal_left_cell_score = diagonal_left_cell + match_bonus

                score = max([from_left_score, from_above_score,
                            diagonal_left_cell_score])
                # take the max
                if score < 0:
                    score = 0
                scoring_array[row][col] = score
                # print(score)
                # make note of which cell was the max in the traceback array
                # using Unicode arrows
                if score == 0:
                    arrow = '-'
                elif score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_cell_score:
                    arrow = up_left_arrow
            # print(arrow)
            # print(traceback_array)
            # print(traceback_array[1][1])

            traceback_array[row][col] = arrow
            # print(traceback_array)
            document['tracebackArray'].textContent = [
                ('\n', traceback_array[i]) for i in range(len(traceback_array))]
            document['scoreArray'].textContent = [
                ('\n', scoring_array[i]) for i in range(len(scoring_array))]

            print(score)

            scoring_array[row][col] = score
            print(scoring_array)
    maxScore = max(max(x) for x in scoring_array)
   
    for i in range(len(scoring_array)):
        for j in range(len(scoring_array)):
            if(scoring_array[i][j] == maxScore):
                print(scoring_array[i][j], i, j)
                rowMaxAlign = i
                colMaxAlign = j
                strAlignment = traceback_local_alignment(
                    traceback_array, sequenceFst, sequenceSnd, rowMaxAlign, colMaxAlign)
                break
 
    document['scoredAlignment'].textContent = '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[0], '\n', traceback_alignment(
        traceback_array, sequenceFst, sequenceSnd)[1], '\n', traceback_alignment(traceback_array, sequenceFst, sequenceSnd)[2]


# ____________________________________________________________________________________________________________________

# ____________________________________________________________________________________________________________________
# ____________________________________________________________________________________________________________________


document['global-btn'].bind('click', globalAlignment)
document['globalProtein-btn'].bind('click', globalAlignmentProtein)
document['Autism-ProAlignment-btn'].bind('click', globalAutismAlignmentProtein )
document['Autism-DNA-Alignment-btn'].bind('click', globalAutismAlignmentDNA)
document['localProtein-btn'].bind('click', localProteinAlignment)
document['local-btn'].bind('click', localAlignment)
