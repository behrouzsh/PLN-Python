import json
from search import Client
import sys
reload(sys)
import csv
sys.setdefaultencoding('utf8')

class Pln(object):
    """Pln class is to read the input file and extract the motifs and modifications.
    :param motifs: List of motifs to be searched.
    """
    motif_and_modification_list = []
    search_client = Client()

    def __init__(self):
        self.motif_and_modification_list = []

    def read_config_json_file(self):
        with open('config/config.json', "r") as json_data_file:
            self._config_data = json.load(json_data_file)
        return self._config_data

    def read_input_file(self):
        with open('input/motif.txt', "r") as input_file:
            self._motifs_data = input_file.read().splitlines()
            #print self._motifs_data
        input_file.close()
        return self._motifs_data

    def extract_motifs(self, motifs_data):
        for motif in motifs_data:
            tmp_list = self.analyze_motif(motif)
            self.motif_and_modification_list.append(tmp_list)
        print "----------Here----------"
        print self.motif_and_modification_list
        return self.motif_and_modification_list

    def analyze_motif(self, motif):
        mod_idx = False
        local_motif = ""
        local_modification = ""
        local_diff = ""
        modification_list = []
        mod_place = 0
        idx = 0
        motif_idx = 0
        while idx < len(motif):
            if motif[idx] == "[":
                local_modification = local_modification + motif[idx - 1]
                mod_place = motif_idx
                mod_idx = True

            if mod_idx == False:
                local_motif = local_motif + motif[idx]
                motif_idx = motif_idx + 1

            while mod_idx:
                local_diff = local_diff + motif[idx + 1]
                if motif[idx + 1] == "]":
                    local_diff = local_diff[:-1]
                    mod_idx = False
                    modification_list.append([mod_place,local_modification,local_diff])
                    local_modification = ""
                    local_diff = ""

                idx = idx + 1
            idx = idx + 1
        # print "local_motif =========================="
        # print local_motif
        # print "modification_list =========================="
        # print modification_list
        return [local_motif,modification_list]


    def call_prosite(self):
        self.prosite_result = []
        for item in self.motif_and_modification_list:
            self.prosite_result.append(self.search_client.search_prosite(item[0]))
        return

    def call_psimod(self):
        self.psimod_result = []
        #local_psimod = []
        #local_similar_psimod = []
        item_num = -1
        for item in self.motif_and_modification_list:
            item_num = item_num + 1
            psimod_num = -1
            for psimod in item[1]:
                psimod_num = psimod_num + 1
                with open('pln/resources/psi-mod/mapping.csv') as csvfile:
                    reader = csv.reader(csvfile, delimiter='\t')
                    diff = 100000000
                    first_line = True
                    print "psimod"
                    print psimod[1], '--', psimod[2]  # , '--',row[2], '--',row[3]
                    local_psimod = []
                    local_similar_psimod = []
                    similar = 0
                    for row in reader:
                        if first_line:
                            first_line = False
                            continue
                        #print "row"
                        #print row[0], '--',row[1], '--',row[2]

                        local_diff = self.extract_diff(row[2],psimod[2])
                        #print str(psimod[1]), str(row[1]), str(psimod[1]) == str(row[1])
                        if str(psimod[1]) == str(row[1]) and local_diff < diff:
                            diff = local_diff
                            local_psimod = []
                            local_similar_psimod = []
                            local_psimod.append(row)

                            similar = 0
                            continue
                        if psimod[1] == row[1] and local_diff == diff:
                            local_similar_psimod.append(row)
                            similar = similar + 1
                    print "local_psimod: ", local_psimod
                    print "local_psimod: ", local_psimod[0][0]
                    print self.motif_and_modification_list
                    print self.motif_and_modification_list[item_num][1][psimod_num]
                    self.motif_and_modification_list[item_num][1][psimod_num].append(local_psimod[0][0])
                    print self.motif_and_modification_list[item_num][1][psimod_num]

                    print "Here 2"
                    print self.motif_and_modification_list
                    print "-------------------------------------------"
                    print "Similar: ", similar
                    print "-------------------------------------------"
                    print "local_similar_psimod: ", local_similar_psimod
                    print "++++++++++++++++++++++++++++++++++++++++++++++++++"
                    self.psimod_result.append([local_psimod,similar,local_similar_psimod, psimod[1], psimod[2]])

        print self.psimod_result
        return self.psimod_result

    def extract_diff(self, idx1,idx2):
        return abs(float(idx1) - float(idx2))

    def print_to_file(self, motif_list, prosite_response, psimod_response):
        with open('output/motif_output.txt', "w") as output_file:
            for x in range(0, len(prosite_response)):
                output_file.write(self._motifs_data[x])
                output_file.write(" -> ")
                output_file.write(motif_list[x][0])
                output_file.write(" : ")
                #output_file.write(str(prosite_response[x]))
                output_inchlike = "PLN=ver1:InChl-like/r="
                n_match = prosite_response[x].get('n_match')
                matches = prosite_response[x].get('matchset')

                # uniprot ++++++++++++++++++++++++++++++++++++++++++++
                flag = 0
                for match in range(0, n_match):
                    print matches[match]
                    print "+++++++++"
                    output_inchlike = output_inchlike + "uniprot:" + matches[match].get("sequence_ac")+ ";"
                    flag = 1
                if flag == 1:
                    output_inchlike = output_inchlike[:-1]

                # hugo +++++++++++++++++++++++++++++++++++++++++++++++
                flag = 0
                output_inchlike = output_inchlike + "/s="
                for match in range(0, n_match):
                    output_inchlike = output_inchlike + "hugo:" + matches[match].get("sequence_id") + ";"
                    flag = 1
                if flag == 1:
                    output_inchlike = output_inchlike[:-1]

                # MOD +++++++++++++++++++++++++++++++++++++++++++++++
                output_inchlike = output_inchlike + "/d=/v=/m="
                flag = 0
                for match in range(0, n_match):
                    for modification in range(0, len(self.motif_and_modification_list[x][1])):
                        flag = 1
                        print modification
                        start = matches[match].get("start")
                        print start
                        print int(self.motif_and_modification_list[x][1][modification][0])
                        position = start -1 + int(self.motif_and_modification_list[x][1][modification][0])
                        print position
                        print self.motif_and_modification_list[x][1][modification][3]
                        output_inchlike = output_inchlike + self.motif_and_modification_list[x][1][modification][3] + \
                                          "@" +  str(position) + ";"
                if flag == 1:
                    output_inchlike = output_inchlike[:-1]

                output_inchlike = output_inchlike + "#"

                output_file.write("\n------------------\n")
                output_file.write(output_inchlike)
                output_file.write("\n++++++++++++++++++\n")
            for y in range(0, len(psimod_response)):
                output_file.write("\nModification : %s[%s]\n" %(str(psimod_response[y][3]), str(psimod_response[y][4])))
                output_file.write(str(psimod_response[y][0]))
                output_file.write("\nNumber of similars : \n")
                output_file.write(str(psimod_response[y][1]))
                output_file.write("\nSimilar : \n")
                for similar in psimod_response[y][2]:
                    output_file.write(str(similar))
                    output_file.write("\n")
                output_file.write("\n-------------------\n")

        output_file.close()

        return




