
## Importing packages
import numpy as np
import pandas as pd
import xlsxwriter

## Reading the relevant sysmex excel columns into a Pandas datafile
file_loc = r'C:\Users\mogmi\OneDrive\Documents\Medisin\Prosjektoppgave\sysmex_anonymous.xlsx'
columns_for_use = 'B, J, BM, BQ, BV, BW, BX, CB, CF, CL, CP, DB, DD, DF, DH, DJ, DL, DV, DZ, ED'
df = pd.read_excel(file_loc, index_col=None, na_values=['NA'], usecols = columns_for_use)


## Naming the columns
df.columns = [
    'days', 'patient_numbers', 'flagblast', 'flagatypical', 'flagfrag',
    'flagclump', 'WBC', 'HGB', 'MCV', 'PLT',
    'RDW', 'NRBC', 'neutro', 'lympho', 'mono',
    'eos', 'baso', 'basopros', 'immatureg', 'retic']
days = df['days']

## Making a list of the patient numbers without duplicates
patient_numbers_excel = df['patient_numbers']
patient_numbers = list(dict.fromkeys(patient_numbers_excel)) # patient numbers without duplicates
n_patients = len(patient_numbers)

print('Number of patients:')
print(n_patients)
np.savetxt('pasientnummere.txt', patient_numbers, '%6d')



## Makingn matrix for tracking each patient's test location in the excel file
first_and_previous = np.zeros( (n_patients, 85) ) # Zero-matrix meant for the excel location of the first test after 48 days, and the previous
first_and_previous[:,0] = patient_numbers # Sets column 0 as patient number column


## Finding first test after 48 days, plus previous tests for same patient.
cutoff = 48 # cutoff for days
for i in range( len(days) ): # i is the same as the excel index
    if (i > 0) and (days[i] >= cutoff): # If "day" >= cutoff
        if (days[i-1] < cutoff) and (patient_numbers_excel[i-1] == patient_numbers_excel[i]): # If previous "day" < cutoff AND it's the same patient (= if patient's first day over cutoff, but has previous days)
            patientnumber_index = patient_numbers.index(patient_numbers_excel[i]) #Finding the index of the patient in the patient number list not containing duplicates
            j = 0
            while (patient_numbers_excel[i-j] == patient_numbers_excel[i]): # While it's the same patient
                first_and_previous[patientnumber_index, j+1] = i - j # Assigning the excel x-axil index of the test, and previous to the "irst and prev" (+2 would make excel indexing)
                j += 1
        if (patient_numbers_excel[i-1] != patient_numbers_excel[i]): # If previous is not the same patient
            patientnumber_index = patient_numbers.index(patient_numbers_excel[i])
            first_and_previous[patientnumber_index, 1] = i # (+2 would make excel indexing)

## Creating a more compact matrix, excluding patients without tests after the cutoff.
compact_firstandprev = first_and_previous[first_and_previous[:,1] != 0]
compact_firstandprev2 = compact_firstandprev[:,compact_firstandprev.any(axis=0)]
print('Number of patients with tests after cutoff (compact_firstandprev2):')
print(np.shape(compact_firstandprev2))
np.savetxt('tests_file.txt', compact_firstandprev, '%5d')



### LRGR ### - Algorithm for applying all "> -rules"

lrgr_matrix = np.zeros( (4295,20) ) # 4295 is the length of compact_firstandprev 1 and 2, and 5 is for flags: RBC, PLT, atypisk, Blast
lrgr_matrix[:,0] = compact_firstandprev[:,0]

          #       1      2      3      4      5       6         7         8      9      10         11      12      13      14          15          16          17             18            19
input_column = ['WBC', 'PLT', 'HGB', 'MCV', 'RDW', 'neutro', 'lympho', 'mono', 'eos', 'baso', 'basopros', 'NRBC', 'NRBC', 'retic', 'immatureg', 'flagfrag', 'flagclump', 'flagatypical', 'flagblast']
input_limit =  [30,    1000,   19,    105,   22,      20,       5,      1.5,     2,    0.5,       3,        1,     0.001,   100,       5,          99,         99,         99,           99]

for n in range(len(input_column)): # iterates over input column and its limits. (19)
    for i in range(4295): # iterates over patients
        stopp = 0 # Stopp means that flag/output is found
        test_index = compact_firstandprev[i,1]
        element = df.at[test_index, input_column[n]]
        if isinstance(element,str) != True: # moves on if not string
            if element > input_limit[n]: # If this patient's first test is flagged
                lrgr_matrix[i,n+1] = 2 # 2 = "flag"
                if compact_firstandprev[i,2] == 0: # If there are no previous tests
                    stopp = 1
            if element < input_limit[n]: #Stops if no flag is found
                stopp = 1
            j = 2 # Starting by looking at previous test
            while stopp == 0 and compact_firstandprev[i,j] != 0:
                prev_test_index = compact_firstandprev[i,j]
                element_prev = df.at[prev_test_index, input_column[n]]
                if isinstance(element_prev,str) != True: #Moves on if not string
                    if element_prev > input_limit[n]:
                        lrgr_matrix[i,n+1] = 1 # 1 = previous has flags
                        stopp = 1
                j += 1

np.savetxt('lrgr_file.txt', lrgr_matrix, '%5d')


### SMLR ### - Algorithm for applying all "< -rules"

smlr_matrix = np.zeros( (4295,7) ) # 4295 is the length of compact_firstandprev 1 and 2, and 5 is for flags: RBC, PLT, atypisk, Blast
smlr_matrix[:,0] = compact_firstandprev[:,0]

          #       1      2      3      4      5       6
input_column = ['WBC', 'WBC', 'PLT', 'HGB', 'MCV', 'neutro']
input_limit =   [3,      4,    100,    7,    75,      1]

for n in range(len(input_column)): # iterates over input column and its limits. (19)
    for i in range(4295): # iterates over patients
        stopp = 0 # Stopp means that flag/output is found
        test_index = compact_firstandprev[i,1]
        element = df.at[test_index, input_column[n]]
        if isinstance(element,str) != True: # moves on if not string
            if element < input_limit[n]: # If this patient's first test is flagged
                smlr_matrix[i,n+1] = 2 # 2 = "flag"
                if compact_firstandprev[i,2] == 0: # If there are no previous tests
                    stopp = 1
            if element > input_limit[n]: # Stops if no flag is found
                stopp = 1
            j = 2 # Starting by looking at previous test
            while stopp == 0 and compact_firstandprev[i,j] != 0:
                prev_test_index = compact_firstandprev[i,j]
                element_prev = df.at[prev_test_index, input_column[n]]
                if isinstance(element_prev,str) != True: # Moves on if not string
                    if element_prev < input_limit[n]:
                        smlr_matrix[i,n+1] = 1 # 1 = previous has flags
                        stopp = 1
                j += 1


np.savetxt('smlr_file.txt', smlr_matrix, '%5d')



## Making matrices for recording rule satisfaction
rule_satisfaction = np.zeros( (4295,30) ) # 17 + 13 = 30 rules
rule_satisfaction[:,0] = compact_firstandprev[:,0]

rule_statistics = np.zeros( (2,29) )
rule_statistics[0,:] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 5, 7, 9, 10, 15, 17, 18, 19, 20, 21, 22, 23, 27, 30, 35, 36]
           # rule_n =   1  2  3  4  5  6  7  8  9  10  11  12  13 14 15 16  17  18  19  20  21  22  23  24  25  26  27  28  29



### IMPLEMENTING RULES ###

## Lrgr (1 og 2)
rules_lrgr =    [3,     11,          12,          13,    26,       27]
              # RDW, RBC_frag, atypiske_lymfo, blast, RBC_frag, PLT_agg
lrgr_columns_1 = [5,    16,          18,          19,    16,       17]
for i in range(4295): # Iterates over patients numbers (row)
    for j in range(len(rules_lrgr)): # Iterates over
        if lrgr_matrix[ i, lrgr_columns_1[j] ] >= 1:
            rule_satisfaction[i, rules_lrgr[j]] = 1

## Lrgr first time (2)
rules_lrgr_ft = [1,   2,   5,     6,    7,   9,   10,     14,   15, 16,  17,  18,    19,     20,   21,  22,   23,   24,   25,       28,         29]
              # WBC, PLT, lymfo, mono, eos, NRBC, umodne, WBC, PLT, HGB, MCV, RDW, neutro, lymfo, mono, eos, baso, NRBC, retic, atypisk_lymfo, blast
lrgr_columns_2 = [1,  2,   7,     8,    9,   12,  15,     1,    2,   3,   4,   5,    6,      7,    8,    9,   10,   13,   14,       18,         19]
for i in range(4295):
    for j in range(len(rules_lrgr_ft)):
        if lrgr_matrix[ i, lrgr_columns_2[j] ] == 2:
            rule_satisfaction[i, rules_lrgr_ft[j]] = 1

## Dobbelkriterie
rules_dobbelkriterie = [8]
                # Baso, Basopros
lrgr_columns_3 = [10,     11]
for i in range(4295):
    if lrgr_matrix[ i, 10 ] == 2 or lrgr_matrix[ i, 11 ] == 2:
        rule_satisfaction[ i, 8 ] = 1

## Smlr first time
rules_smlr_ft = [1,   2,    4,    14,  15,  16,  17,   19]
              # WBC, PLT, neutro, WBC, PLT, HGB, MCV, neutro
smlr_columns = [1,    2,    6,     1,  2,    3,   4,    6]
for i in range(4295):
    for j in range(len(rules_smlr_ft)):
        if lrgr_matrix[ i, smlr_columns[j] ] == 2:
            rule_satisfaction[i, rules_smlr_ft[j]] = 1



for i in range(29):
    rule_statistics[1,i] = sum(rule_satisfaction[:,i+1])
np.savetxt('rule_statistics.txt', rule_statistics, '%5d')



pos_tests = np.zeros( (4295,2) )
pos_tests[:,0] = compact_firstandprev[:,0]
for i in range(4295):
    if sum(rule_satisfaction[i,1:]) > 0:
        pos_tests[i,1] = 1


print("Positive tests:")
print(sum(pos_tests[:,1]))

print("Shape of rule_satisfaction:")
print(np.shape(rule_satisfaction))





workbook = xlsxwriter.Workbook('pos_tests_table.xlsx')
worksheet = workbook.add_worksheet()

array = rule_satisfaction

row = 0

for col, data in enumerate(array):
    worksheet.write_column(row, col, data)

workbook.close()
