from random import randint
l = 10
ADN = "ATGC"
fasta = open("fastag1.fq","w")
fasta.write(">1"+'\n')
for i in range(0,l):
    hasard = randint(0,3)
    fasta.write(ADN[hasard])
fasta.close()
fasta = open("fastag2.fq","w")
fasta.write(">1"+'\n')
for i in range(0,l):
    hasard = randint(0,3)
    fasta.write(ADN[hasard])
fasta.close()
