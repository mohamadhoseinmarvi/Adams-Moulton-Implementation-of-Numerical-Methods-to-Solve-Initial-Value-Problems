import methods

def main():
    i = open('input.txt', 'r')
    o = open('output.txt', 'w')
    o.close()
    for line in i:
        cmd = line.split()
        method = getattr(methods, cmd[0])
        method(cmd[1:len(cmd)])
    i.close()

if __name__ == '__main__':
    main()