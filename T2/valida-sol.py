#!/usr/bin/python
# coding=utf8

import sys

def carregaInstancia(arq, erro):
# le primeira linha
  linha=arq.readline()
  pedacos=linha.split(' ')
  n=int(pedacos[0])
  m=int(pedacos[1])
  C=int(pedacos[2])
  print("n=%d m=%d C=%d" % (n,m,C))
# le valores dos itens
  linha=arq.readline()
  pedacos=linha.strip().split(' ')
  print(pedacos)
  valor = [ int(x) for x in pedacos ]
  for i in range(0,n):
    print("valor[%d]=%d" % (i+1,valor[i]))
# le pesos dos elementos
  linha=arq.readline()
  pedacos=linha.strip().split(' ')
  peso = [ int(x) for x in pedacos ]
  for i in range(0,m):
    print("peso[%d]=%d" % (i+1,peso[i]))
# le relacao R[i][j]
  R = [[]] * n
  for i in range(0,n):
    linha=arq.readline()
    pedacos=linha.strip().split(' ')
    R[i] = [ int(x) for x in pedacos ]
  for i in range(0,n):
    print("\nconjunto do item %d :" % (i+1), end="")
    for j in range(0,m):
      if R[i][j]>0:
        print("%d " % (j+1), end="")
  return (n,m,C,valor,peso,R)
def pertence(k,conjunto):
  for i in range(0,len(conjunto)):
    if k==conjunto[i]:
      return 1
  return 0
# valida a solucao dada no arquivo arq
def valida(arq, n, m, C, valor, peso, R, erro):
# le primeira linha
  linha=arq.readline()
  pedacos=linha.split(' ')
  z=int(pedacos[0])
  somaPeso=int(pedacos[1])
  print("\nz=%d somaPeso=%d" % (z,somaPeso))
# le itens
  linha=arq.readline()
  pedacos=linha.strip().split(' ')
  print(pedacos)
  itens = [ int(x) for x in pedacos ]
# valida z
  somaValor=0
  for i in range(0,len(itens)):
    print("%d " % itens[i], end="")
    somaValor+=valor[itens[i]-1]
  if somaValor>z or somaValor<z:
    print("soma dos valores=%d difere de z=%d" % (somaValor,z))
    return "incorreto;soma dos valores dos itens incorreta. Deveria ser %d" % somaValor
# le elementos
  linha=arq.readline()
  pedacos=linha.strip().split(' ')
  print(pedacos)
  elementos = [ int(x) for x in pedacos ]
# marca elementos cobertos pelos itens
  coberto = [0] * m
  for i in range(0,len(itens)):
    for j in range(0,m):
      item = itens[i] - 1
      if R[item][j] and not coberto[j]:
        coberto[j]=1
        if not pertence(j+1,elementos):
          print("Elemento %d foi coberto pelos itens mas nao consta na solucao" % (j+1))
          return "Incorreto;Elemento %d foi coberto pelos itens mas nao consta na solucao" % (j+1)
# valida elementos
  for i in range(0,len(elementos)):
    e=elementos[i]
    if not coberto[e-1]:
      print("elemento %d nao foi coberto pelos itens!" % e)
      return "Incorreto; elemento %d nao foi coberto pelos itens." % e
#valida capacidade
  soma=0
  for i in range(0,len(elementos)):
    print("%d " % elementos[i], end="")
    soma+=peso[elementos[i]-1]
  if soma>somaPeso or soma<somaPeso:
    print("soma dos pesos=%d difere do valor lido=%d" % (soma,somaPeso))
    return("incorreto;soma dos pesos incorreta. Deveria ser %d" % (soma))
  if soma>C:
    return("Soma dos pesos=%d viola a capacidade C=%d" % (soma,C))
  return "correto;"
if len(sys.argv) != 4:
  print("Uso: %s instancia solucao saida")
else:
  arqInst = open(sys.argv[1], "r")
  arqSol = open(sys.argv[2], "r")
  arqOut = open(sys.argv[3], "a")

  (n,m,C,valor,peso,R) = carregaInstancia(arqInst, sys.stderr)
  resposta=valida(arqSol,n,m,C,valor,peso,R,sys.stderr)
  arqOut.write("%s;%s\n" % (sys.argv[2],resposta))
  print("\n%s" % resposta)
