from random import randint
l = 10
AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
fasta = open("fastap1.fq","w")
fasta.write(">1"+'\n')
for i in range(0,l):
    hasard = randint(0,19)
    fasta.write(AA[hasard])
fasta.close()
fasta = open("fastap2.fq","w")
fasta.write(">1"+'\n')
for i in range(0,l):
    hasard = randint(0,19)
    fasta.write(AA[hasard])
fasta.close()
