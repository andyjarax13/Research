def main():
    txt_file = open("clusters.txt").readlines()
    while True:
        if txt_file[0][:4] != "No.|":
            txt_file.pop(0)
        else:
            break
    for i in range(len(txt_file)):
        txt_file[i] = txt_file[i].strip().split("|")
        txt_file[i].pop(0)
        for j in range(len(txt_file[i])):
            txt_file[i][j] = txt_file[i][j].strip()
        txt_file[i] = ",".join(txt_file[i]) + "\n"
    open("clusters.csv", "w").writelines(txt_file)

main()