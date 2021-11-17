count = 0
with open("../../data/Lolium2019_10X_20PHRED_filtered.sync") as f:
    line = f.readline()
    while(line):
        line = line.split()
        line = line[3:]

        include_pop = True
        for population in line:
            population = population.split(":")
            population = [int(e) for e in population]
            #print(sum(population), "", end="")
            if sum(population) < 5:
                include_pop = False
                break
        #print("")
        count += include_pop
        line = f.readline()
print(count)
