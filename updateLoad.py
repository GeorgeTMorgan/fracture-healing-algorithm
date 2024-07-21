def updateLoad(filename, load):

    updatedLoad = '{:.15e}'.format(load)
    updatedLoad = updatedLoad[:-4] + updatedLoad[-3] + updatedLoad[-1]

    with open(filename + '.dat', 'r') as f:
        searchlines = f.readlines()

    for i, line in enumerate(searchlines):
        if "point load" in line: 

            # edit the relevant line in the file
            editLine = searchlines[i+3].split()
            editLine[2] = updatedLoad
            editLine = ' ' + ' '.join(editLine) + '\n'

            # replace the line in the file
            searchlines[i+3] = editLine

    with open(filename + '.dat', 'w') as f:
        updatedFile = ''.join(searchlines)
        f.write(updatedFile)

