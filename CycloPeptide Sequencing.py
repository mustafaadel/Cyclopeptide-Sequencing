from collections import Counter
# Reading the Amino Acids and weights from a file and storing them in a dictionary
def Initial_List(inputspectrum):

  f = open("weight.txt")
  filelist = []
  dict = {}

# passing the weights file into a dictionary
  for line in f:
      line = line.rstrip()
      words = line.split(" ")
      filelist.append(int(words[1]))
      weight = words[1]
      aminoacid = words[0]
      dict[weight] = aminoacid

# comparing the inputspectrum to the weights file
# retrieving the list of 1-mers
  length = 0
  newlist=[]
  initialList = []
  for k in filelist:
      for i in inputspectrum:
          if i == k:
            newlist.append(k)
            length += 1

# Removing Duplicates from 1-mers
  final_list = []
  [final_list.append(x) for x in newlist if x not in final_list]

# creating the list after removal
  for l in dict.keys():
      for k in final_list:
          if int(l)==k:
              initialList.append(dict[l])

  return initialList

# This function checks if the peptide sequence is consistent or not
def isConsistent(subpeptide, spectrumDic):
    linearSpectra = linear_spectrum(subpeptide)
    subpeptidespectrum = Counter(linearSpectra)
    for key in subpeptidespectrum:
        if key not in spectrumDic:
            return False
        if subpeptidespectrum[key] > spectrumDic[key]:
            return False
    return True

# This function extends the consistent temp list with the initail list
def getCombinations(TheMonomers, TheConsistentList):
    Thepermutations = []
    for i in TheConsistentList:
        for j in TheMonomers:
            mer = i + j
            Thepermutations.append(mer)
    return Thepermutations

# This function calculates the weight of the given peptide sequence and returns a list with the AA weights
def linear_spectrum(peptide):
    f = open("weight.txt")
    dictionaryofAa = {}
    listtofAa = []
    for line in f:
        line = line.rstrip()
        words = line.split(" ")
        listtofAa.append(int(words[1]))
        weights = words[1]
        aminoacid = words[0]
        dictionaryofAa[aminoacid] = weights
    subPeptideList = []
    for i in range(len(peptide)):
        for j in range(i, len(peptide)):
            subPeptideList.append(peptide[i:j + 1])
    subPeptideWeights = []
    for i in subPeptideList:
        total = 0
        for j in i:
            x = int(dictionaryofAa[j])
            total = total + x
        subPeptideWeights.append(total)
    subPeptideWeights.sort()
    print(subPeptideList)
    return subPeptideWeights

# Main
#theoriticalSpectrum = list(input("please input the spectrum list: ").split(" "))
# Convert string to int using list comprehension
#theoriticalSpectrum = [int(i) for i in theoriticalSpectrum]
theoriticalSpectrum = [0 ,97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295, 297, 299, 299, 301, 394, 396, 398, 400, 400, 497]
Initial_L = Initial_List(theoriticalSpectrum)

spectrumDic = Counter(theoriticalSpectrum)

# Count frequency of every char
InitialDic = Counter(Initial_L)

TempList = Initial_L.copy()

# Stop looping when tempList is empty as it carries all the consistant meres
for m in range(len(Initial_L)):
    TempList = getCombinations(Initial_L, TempList)
    consistantList = []
    for i in TempList:
        if isConsistent(i,spectrumDic):
            consistantList.append(i)
    # Stop when all are inconsistant
    if len(consistantList) == 0:
        break
    TempList = set(consistantList.copy())

finalList = TempList
print(finalList)
