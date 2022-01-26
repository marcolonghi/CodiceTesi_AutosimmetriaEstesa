import os
import re
    
yourpath = 'pathPla/'

for root, dirs, files in os.walk(yourpath, topdown=False):
    for name in files:
        print(os.path.join(name))
        nome = os.path.join(name)
        nome = nome.replace(".pla", "")

        with open(yourpath + os.path.join(name), "r") as file: 
            lines = file.readlines()

        for line in lines:
            if(line[0]=="." and line[1]=="i"):
                inp = re.findall(".i\s(\d+)", line)
                inp = str(inp[0])
                inp = int(inp)
            elif(line[0]=="." and line[1]=="o"):
                out = re.findall(".o\s(\d+)", line)
                out = str(out[0])
                out = int(out)
                testo = ".i " + str(inp) + "\n" + ".o 1" + "\n"
                for i in range(1, out+1):
                    f = open("plaCirianiS_split/" + nome + "_" + str(i) + ".pla", "w")
                    f.write(testo)
                    f.close()
            elif(line[0]=="." and line[1]=="e"):
                print("Finito")
                for i in range(1, out+1):
                    f = open("plaCirianiS_split/" + nome + "_" + str(i) + ".pla", "a")
                    f.write(".e")
                    f.close()
            else:
                mint = re.findall("([0,1,-]+)\s[0,1,-]+", line)
                mint = str(mint[0])
                mintout = re.findall("[0,1,-]+\s([0,1,-]+)", line)
                for i in range(0, out):
                    if(mintout[0][i] == "1"):
                        f = open("plaCirianiS_split/" + nome + "_" + str(i+1) + ".pla", "a")
                        f.write(mint + " 1" + "\n")
                        f.close()
                    elif(mintout[0][i] == "-"):
                        f = open("plaCirianiS_split/" + nome + "_" + str(i+1) + ".pla", "a")
                        f.write(mint + " -" + "\n")
                        f.close()
