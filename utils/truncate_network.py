import constants
import networkx as nx
import os
import pandas as pd
import os

import networkx as nx
import pandas as pd

import constants


def read_network(network_file_name):
    edges=[]
    df_network=pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_file_name), sep='\t')
    for i, cur in df_network.iterrows():
        # if i==10000: break
        edges.append((cur[0], cur[2]))

    G = nx.Graph()
    G.add_edges_from(edges)

    return G, df_network.columns

def truncate_network_by_gene_list(gene_list, network_file_name, radius=5):
    nw, headers=read_network(network_file_name)

    total_edges=[]
    for cur_gene in gene_list:
        if cur_gene in nw.nodes():
            total_edges+=list(nx.ego_graph(nw, cur_gene, radius=radius).edges())


    new_network=[]
    for cur_edge in total_edges:
        srt=sorted([cur_edge[0], cur_edge[1]])
        new_network.append("{}\t{}\t{}".format(srt[0], "ppi", srt[1]))

    network_output="\t".join(headers) + "\n"
    network_output+="\n".join(list(set(new_network)))
    file(os.path.join(constants.NETWORKS_DIR, network_file_name[:-4]+"_truncated_{}.sif".format(radius)),'w+').write(network_output)


if __name__=="__main__":
    gene_list = ["ENSG00000146674" ,"ENSG00000105825" ,"ENSG00000100345" ,"ENSG00000107796" ,"ENSG00000072110" ,"ENSG00000142871" ,"ENSG00000196924" ,"ENSG00000163017" ,"ENSG00000108821" ,"ENSG00000187514" ,"ENSG00000130402" ,"ENSG00000106366" ,"ENSG00000075624" ,"ENSG00000210082" ,"ENSG00000149591" ,"ENSG00000166923" ,"ENSG00000182718" ,"ENSG00000122861" ,"ENSG00000167460" ,"ENSG00000103257" ,"ENSG00000128591" ,"ENSG00000095752" ,"ENSG00000117289" ,"ENSG00000198467" ,"ENSG00000167900" ,"ENSG00000176619" ,"ENSG00000187498" ,"ENSG00000135480" ,"ENSG00000128564" ,"ENSG00000205542" ,"ENSG00000131473" ,"ENSG00000197451" ,"ENSG00000101335" ,"ENSG00000117399" ,"ENSG00000210196" ,"ENSG00000166508" ,"ENSG00000184009" ,"ENSG00000159167" ,"ENSG00000177469" ,"ENSG00000176014" ,"ENSG00000123416" ,"ENSG00000091136" ,"ENSG00000092820" ,"ENSG00000071127" ,"ENSG00000231991" ,"ENSG00000120129" ,"ENSG00000198695" ,"ENSG00000196419" ,"ENSG00000167553" ,"ENSG00000187837" ,"ENSG00000130513" ,"ENSG00000180573" ,"ENSG00000140416" ,"ENSG00000164111" ,"ENSG00000135069" ,"ENSG00000249992" ,"ENSG00000169429" ,"ENSG00000092621" ,"ENSG00000130176" ,"ENSG00000167157" ,"ENSG00000004799" ,"ENSG00000110880" ,"ENSG00000149257" ,"ENSG00000178919" ,"ENSG00000213626" ,"ENSG00000064666" ,"ENSG00000165280" ,"ENSG00000164442" ,"ENSG00000204262" ,"ENSG00000102265" ,"ENSG00000108518" ,"ENSG00000106211" ,"ENSG00000125844" ,"ENSG00000197457" ,"ENSG00000181019" ,"ENSG00000163346" ,"ENSG00000176046" ,"ENSG00000083444" ,"ENSG00000108106" ,"ENSG00000197746" ,"ENSG00000002834" ,"ENSG00000103187" ,"ENSG00000184678" ,"ENSG00000152952" ,"ENSG00000198830" ,"ENSG00000152661" ,"ENSG00000136026" ,"ENSG00000188486" ,"ENSG00000160789" ,"ENSG00000188229" ,"ENSG00000176170" ,"ENSG00000123975" ,"ENSG00000198363" ,"ENSG00000059804" ,"ENSG00000148926" ,"ENSG00000034063" ,"ENSG00000137076" ,"ENSG00000159176" ,"ENSG00000088325" ,"ENSG00000176788" ,"ENSG00000116649" ,"ENSG00000143384" ,"ENSG00000101057" ,"ENSG00000147065" ,"ENSG00000183087" ,"ENSG00000134057" ,"ENSG00000159840" ,"ENSG00000110092" ,"ENSG00000112118" ,"ENSG00000136859" ,"ENSG00000096384" ,"ENSG00000110108" ,"ENSG00000162704" ,"ENSG00000187840" ,"ENSG00000112697" ,"ENSG00000181104" ,"ENSG00000026025" ,"ENSG00000198840" ,"ENSG00000123989" ,"ENSG00000162909" ,"ENSG00000157613" ,"ENSG00000105202" ,"ENSG00000260032" ,"ENSG00000101412" ,"ENSG00000130635" ,"ENSG00000124762" ,"ENSG00000154096" ,"ENSG00000164032" ,"ENSG00000067057" ,"ENSG00000109971" ,"ENSG00000085117" ,"ENSG00000099860" ,"ENSG00000163430" ,"ENSG00000168615" ,"ENSG00000134333" ,"ENSG00000057019" ,"ENSG00000171223" ,"ENSG00000060718" ,"ENSG00000153234" ,"ENSG00000169908" ,"ENSG00000143549" ,"ENSG00000150938" ,"ENSG00000150991" ,"ENSG00000164096" ,"ENSG00000166851" ,"ENSG00000073111" ,"ENSG00000117318" ,"ENSG00000011304" ,"ENSG00000171951" ,"ENSG00000090861" ,"ENSG00000089280" ,"ENSG00000130066" ,"ENSG00000205426" ,"ENSG00000131016" ,"ENSG00000132646" ,"ENSG00000137801" ,"ENSG00000125166" ,"ENSG00000065978" ,"ENSG00000103241" ,"ENSG00000065308" ,"ENSG00000197879" ,"ENSG00000172216" ,"ENSG00000168003" ,"ENSG00000235162" ,"ENSG00000183684" ,"ENSG00000116133" ,"ENSG00000128606" ,"ENSG00000254851" ,"ENSG00000135919" ,"ENSG00000099783" ,"ENSG00000188153" ,"ENSG00000120738" ,"ENSG00000158373" ,"ENSG00000197043" ,"ENSG00000128595" ,"ENSG00000198786" ,"ENSG00000143799" ,"ENSG00000122884" ,"ENSG00000100297" ,"ENSG00000068697" ,"ENSG00000108298" ,"ENSG00000137124" ,"ENSG00000172757" ,"ENSG00000184374" ,"ENSG00000168386" ,"ENSG00000110104" ,"ENSG00000197903" ,"ENSG00000131236" ,"ENSG00000182481" ,"ENSG00000072041" ,"ENSG00000170775" ,"ENSG00000128272" ,"ENSG00000164104" ,"ENSG00000035403" ,"ENSG00000133657" ,"ENSG00000130522" ,"ENSG00000171848" ,"ENSG00000196611" ,"ENSG00000107130" ,"ENSG00000187479" ,"ENSG00000105974" ,"ENSG00000165916" ,"ENSG00000115758" ,"ENSG00000256663" ,"ENSG00000179218" ,"ENSG00000087303" ,"ENSG00000159251" ,"ENSG00000183691" ,"ENSG00000132341" ,"ENSG00000074800" ,"ENSG00000142230" ,"ENSG00000129250" ,"ENSG00000127589" ,"ENSG00000196923" ,"ENSG00000167772" ,"ENSG00000132507" ,"ENSG00000169230" ,"ENSG00000168439" ,"ENSG00000091986" ,"ENSG00000124766" ,"ENSG00000168209" ,"ENSG00000125520" ,"ENSG00000112984" ,"ENSG00000076003" ,"ENSG00000100401" ,"ENSG00000163468" ,"ENSG00000135624" ,"ENSG00000108424" ,"ENSG00000169715" ,"ENSG00000133872" ,"ENSG00000067225" ,"ENSG00000120889" ,"ENSG00000117394" ,"ENSG00000111206" ,"ENSG00000125148" ,"ENSG00000162734" ,"ENSG00000050405" ,"ENSG00000125753" ,"ENSG00000163431" ,"ENSG00000167996" ,"ENSG00000131981" ,"ENSG00000168077" ,"ENSG00000145386" ,"ENSG00000175040" ,"ENSG00000198901" ,"ENSG00000143013" ,"ENSG00000204397" ,"ENSG00000095303" ,"ENSG00000164171" ,"ENSG00000163739" ,"ENSG00000233328" ,"ENSG00000175334" ,"ENSG00000104738" ,"ENSG00000130147" ,"ENSG00000104889" ,"ENSG00000105323" ,"ENSG00000139329" ,"ENSG00000103222" ,"ENSG00000134013" ,"ENSG00000117519" ,"ENSG00000117450" ,"ENSG00000128245" ,"ENSG00000178209" ,"ENSG00000129195" ,"ENSG00000104368" ,"ENSG00000078369" ,"ENSG00000175063" ,"ENSG00000111897" ,"ENSG00000181649" ,"ENSG00000113369" ,"ENSG00000107562" ,"ENSG00000123496" ,"ENSG00000138061" ,"ENSG00000134690" ,"ENSG00000168496" ,"ENSG00000104904" ,"ENSG00000067560" ,"ENSG00000013297" ,"ENSG00000101447" ,"ENSG00000141429" ,"ENSG00000184640" ,"ENSG00000105281" ,"ENSG00000197744" ,"ENSG00000116729" ,"ENSG00000122862" ,"ENSG00000130726" ,"ENSG00000087586" ,"ENSG00000013810" ,"ENSG00000130816" ,"ENSG00000167513" ,"ENSG00000260549" ,"ENSG00000092199" ,"ENSG00000113368" ,"ENSG00000100292" ,"ENSG00000114867" ,"ENSG00000155760" ,"ENSG00000198886" ,"ENSG00000173457" ,"ENSG00000147224" ,"ENSG00000099901" ,"ENSG00000136450" ,"ENSG00000198431" ,"ENSG00000177706" ,"ENSG00000106628" ,"ENSG00000010292" ,"ENSG00000073008" ,"ENSG00000100234" ,"ENSG00000104321" ,"ENSG00000109321" ,"ENSG00000143321" ,"ENSG00000166986" ,"ENSG00000109099" ,"ENSG00000149090" ,"ENSG00000090273" ,"ENSG00000004776" ,"ENSG00000156515" ,"ENSG00000159335" ,"ENSG00000076513" ,"ENSG00000184216" ,"ENSG00000138668" ,"ENSG00000102359" ,"ENSG00000115884" ,"ENSG00000129116" ,"ENSG00000070669" ,"ENSG00000176907" ,"ENSG00000178999" ,"ENSG00000101182" ,"ENSG00000175592" ,"ENSG00000206053" ,"ENSG00000166741" ,"ENSG00000143621" ,"ENSG00000136156" ,"ENSG00000171793" ,"ENSG00000092841" ,"ENSG00000154277" ,"ENSG00000113758" ,"ENSG00000111669" ,"ENSG00000124702" ,"ENSG00000134531" ,"ENSG00000160211" ,"ENSG00000254999" ,"ENSG00000170515" ,"ENSG00000244486" ,"ENSG00000180914" ,"ENSG00000125835" ,"ENSG00000116478" ,"ENSG00000122952" ,"ENSG00000168393" ,"ENSG00000133816" ,"ENSG00000158710" ,"ENSG00000013306" ,"ENSG00000166147" ,"ENSG00000075618" ,"ENSG00000115641" ,"ENSG00000074695" ,"ENSG00000110107" ,"ENSG00000196262" ,"ENSG00000214706" ,"ENSG00000149503" ,"ENSG00000140682" ,"ENSG00000176890" ,"ENSG00000175550" ,"ENSG00000173726" ,"ENSG00000063660" ,"ENSG00000094804" ,"ENSG00000115053" ,"ENSG00000019991" ,"ENSG00000198888" ,"ENSG00000196526" ,"ENSG00000188612" ,"ENSG00000129351" ,"ENSG00000117395" ,"ENSG00000136938" ,"ENSG00000115457" ,"ENSG00000007080" ,"ENSG00000138069" ,"ENSG00000131462" ,"ENSG00000181218" ,"ENSG00000126457" ,"ENSG00000179820" ,"ENSG00000165283" ,"ENSG00000157456" ,"ENSG00000104852" ,"ENSG00000149136" ,"ENSG00000102760" ,"ENSG00000136158" ,"ENSG00000087087" ,"ENSG00000198176" ,"ENSG00000108294" ,"ENSG00000079616" ,"ENSG00000106799" ,"ENSG00000102144" ,"ENSG00000075785" ,"ENSG00000006016" ,"ENSG00000137168" ,"ENSG00000134684" ,"ENSG00000138119" ,"ENSG00000119139" ,"ENSG00000101361" ,"ENSG00000105011" ,"ENSG00000104611" ,"ENSG00000143418" ,"ENSG00000087086" ,"ENSG00000163902" ,"ENSG00000174807" ,"ENSG00000110917" ,"ENSG00000130204" ,"ENSG00000135318" ,"ENSG00000088247" ,"ENSG00000071539" ,"ENSG00000143761" ,"ENSG00000116717" ,"ENSG00000127920" ,"ENSG00000020181" ,"ENSG00000025770" ,"ENSG00000143545" ,"ENSG00000105223" ,"ENSG00000089597" ,"ENSG00000161203" ,"ENSG00000108561" ,"ENSG00000048392" ,"ENSG00000108829" ,"ENSG00000077312" ,"ENSG00000211459" ,"ENSG00000173207" ,"ENSG00000123136" ,"ENSG00000127528" ,"ENSG00000187678" ,"ENSG00000161326" ,"ENSG00000140545" ,"ENSG00000058668" ,"ENSG00000142453" ,"ENSG00000165169" ,"ENSG00000093009" ,"ENSG00000206625" ,"ENSG00000108604" ,"ENSG00000106537" ,"ENSG00000116221" ,"ENSG00000129103" ,"ENSG00000170558" ,"ENSG00000150753" ,"ENSG00000198959" ,"ENSG00000240342" ,"ENSG00000135862" ,"ENSG00000116560" ,"ENSG00000142945" ,"ENSG00000142949" ,"ENSG00000183207" ,"ENSG00000005022" ,"ENSG00000167779" ,"ENSG00000106244" ,"ENSG00000033050" ,"ENSG00000164611" ,"ENSG00000075213" ,"ENSG00000100162" ,"ENSG00000062822" ,"ENSG00000136810" ,"ENSG00000119669" ,"ENSG00000077152" ,"ENSG00000128283" ,"ENSG00000166073" ,"ENSG00000162512" ,"ENSG00000187653" ,"ENSG00000166340" ,"ENSG00000100644" ,"ENSG00000166963" ,"ENSG00000170525" ,"ENSG00000125534" ,"ENSG00000003436" ,"ENSG00000110047" ,"ENSG00000119681" ,"ENSG00000044115" ,"ENSG00000215301" ,"ENSG00000148484" ,"ENSG00000004142" ,"ENSG00000187990" ,"ENSG00000229132" ,"ENSG00000167470" ,"ENSG00000234743" ,"ENSG00000060138" ,"ENSG00000161800" ,"ENSG00000120254" ,"ENSG00000136997" ,"ENSG00000097021" ,"ENSG00000155660" ,"ENSG00000004478" ,"ENSG00000137575" ,"ENSG00000130985" ,"ENSG00000136699" ,"ENSG00000121152" ,"ENSG00000130203" ,"ENSG00000143368" ,"ENSG00000005884" ,"ENSG00000185885" ,"ENSG00000126254" ,"ENSG00000167670" ,"ENSG00000146425" ,"ENSG00000198727" ,"ENSG00000106682" ,"ENSG00000179967" ,"ENSG00000197785" ,"ENSG00000181163" ,"ENSG00000139636" ,"ENSG00000063244" ,"ENSG00000115738" ,"ENSG00000025772" ,"ENSG00000110955" ,"ENSG00000075218" ,"ENSG00000143575" ,"ENSG00000116209" ,"ENSG00000166710" ,"ENSG00000241685" ,"ENSG00000183963" ,"ENSG00000198732" ,"ENSG00000182934" ,"ENSG00000100029" ,"ENSG00000159674" ,"ENSG00000167601" ,"ENSG00000175567" ,"ENSG00000163950" ,"ENSG00000175166" ,"ENSG00000118523" ,"ENSG00000141401" ,"ENSG00000183688" ,"ENSG00000139112" ,"ENSG00000140859" ,"ENSG00000139514" ,"ENSG00000125384" ,"ENSG00000168476" ,"ENSG00000068489" ,"ENSG00000126950" ,"ENSG00000049541" ,"ENSG00000169857" ,"ENSG00000134352" ,"ENSG00000132382" ,"ENSG00000124562" ,"ENSG00000265415" ,"ENSG00000161544" ,"ENSG00000024422" ,"ENSG00000154734" ,"ENSG00000145912" ,"ENSG00000130706" ,"ENSG00000134294" ,"ENSG00000076706" ,"ENSG00000133275" ,"ENSG00000114767" ,"ENSG00000169710" ,"ENSG00000141905" ,"ENSG00000185650" ,"ENSG00000126945" ,"ENSG00000102038" ,"ENSG00000112081" ,"ENSG00000165502" ,"ENSG00000031698" ,"ENSG00000130340" ,"ENSG00000127564" ,"ENSG00000053372" ,"ENSG00000171206" ,"ENSG00000168906" ,"ENSG00000182199" ,"ENSG00000140564" ,"ENSG00000163466" ,"ENSG00000133110" ,"ENSG00000105220" ,"ENSG00000116237" ,"ENSG00000198961" ,"ENSG00000108953" ,"ENSG00000074071" ,"ENSG00000140285" ,"ENSG00000114251" ,"ENSG00000142173" ,"ENSG00000115875" ,"ENSG00000069956" ,"ENSG00000161547" ,"ENSG00000153214" ,"ENSG00000180340" ,"ENSG00000160685" ,"ENSG00000115380" ,"ENSG00000179115" ,"ENSG00000146670" ,"ENSG00000102531" ,"ENSG00000130508" ,"ENSG00000034152" ,"ENSG00000119801" ,"ENSG00000187601" ,"ENSG00000076382" ,"ENSG00000174775" ,"ENSG00000115363" ,"ENSG00000206190" ,"ENSG00000225485" ,"ENSG00000163659" ,"ENSG00000112769" ,"ENSG00000111640" ,"ENSG00000151233" ,"ENSG00000196954" ,"ENSG00000111602" ,"ENSG00000102409" ,"ENSG00000179262" ,"ENSG00000203760" ,"ENSG00000203812" ,"ENSG00000128050" ,"ENSG00000077549" ,"ENSG00000183558" ,"ENSG00000198858" ,"ENSG00000177156" ,"ENSG00000074370" ,"ENSG00000188976" ,"ENSG00000185033" ,"ENSG00000106025" ,"ENSG00000135486" ,"ENSG00000185624" ,"ENSG00000079257" ,"ENSG00000256269" ,"ENSG00000166226" ,"ENSG00000182220" ,"ENSG00000101220" ,"ENSG00000198755" ,"ENSG00000103653" ,"ENSG00000168066" ,"ENSG00000179091" ,"ENSG00000157227" ,"ENSG00000136238" ,"ENSG00000136942" ,"ENSG00000100316" ,"ENSG00000141696" ,"ENSG00000134259" ,"ENSG00000213465" ,"ENSG00000121774" ,"ENSG00000128342" ,"ENSG00000119333" ,"ENSG00000196154" ,"ENSG00000100714" ,"ENSG00000079308" ,"ENSG00000111799" ,"ENSG00000137207" ,"ENSG00000138434" ,"ENSG00000109846" ,"ENSG00000179041" ,"ENSG00000115275" ,"ENSG00000115919" ,"ENSG00000100304" ,"ENSG00000104884" ,"ENSG00000166582" ,"ENSG00000168090" ,"ENSG00000111641" ,"ENSG00000116017" ,"ENSG00000162458" ,"ENSG00000164733" ,"ENSG00000131467" ,"ENSG00000106263" ,"ENSG00000172935" ,"ENSG00000170801" ,"ENSG00000166825" ,"ENSG00000112658" ,"ENSG00000124882" ,"ENSG00000147164" ,"ENSG00000184900" ,"ENSG00000142507" ,"ENSG00000180198" ,"ENSG00000086758" ,"ENSG00000104765" ,"ENSG00000101444" ,"ENSG00000085662" ,"ENSG00000253368" ,"ENSG00000172301" ,"ENSG00000138107" ,"ENSG00000110002" ,"ENSG00000106105" ,"ENSG00000065427" ,"ENSG00000167641" ,"ENSG00000070404" ,"ENSG00000108679" ,"ENSG00000011009" ,"ENSG00000055070" ,"ENSG00000170144" ,"ENSG00000139318" ,"ENSG00000228716" ,"ENSG00000197565" ,"ENSG00000158406" ,"ENSG00000144354" ,"ENSG00000240476" ,"ENSG00000109272" ,"ENSG00000139641" ,"ENSG00000168028" ,"ENSG00000108946" ,"ENSG00000135111" ,"ENSG00000134222" ,"ENSG00000221743" ,"ENSG00000101210" ,"ENSG00000137845" ,"ENSG00000182944" ,"ENSG00000160877" ,"ENSG00000176946" ,"ENSG00000212907" ,"ENSG00000077147" ,"ENSG00000103335" ,"ENSG00000138675" ,"ENSG00000136888" ,"ENSG00000174231" ,"ENSG00000171724" ,"ENSG00000013275" ,"ENSG00000118705" ,"ENSG00000066044" ,"ENSG00000141424" ,"ENSG00000102119" ,"ENSG00000126803" ,"ENSG00000164938" ,"ENSG00000213763" ,"ENSG00000128951" ,"ENSG00000116871" ,"ENSG00000105968" ,"ENSG00000114353" ,"ENSG00000111328" ,"ENSG00000138363" ,"ENSG00000110422" ,"ENSG00000165030" ,"ENSG00000125970" ,"ENSG00000225614" ,"ENSG00000119865" ,"ENSG00000104081" ,"ENSG00000106462" ,"ENSG00000148175" ,"ENSG00000154736" ,"ENSG00000226935" ,"ENSG00000108861" ,"ENSG00000160613" ,"ENSG00000171992" ,"ENSG00000105137" ,"ENSG00000100220" ,"ENSG00000037241" ,"ENSG00000111667" ,"ENSG00000160193" ,"ENSG00000131153" ,"ENSG00000135535" ,"ENSG00000183856" ,"ENSG00000166197" ,"ENSG00000123179" ,"ENSG00000100991" ,"ENSG00000105472" ,"ENSG00000115677" ,"ENSG00000186185" ,"ENSG00000164176" ,"ENSG00000167600" ,"ENSG00000165304" ,"ENSG00000071575" ,"ENSG00000104722" ,"ENSG00000156273" ,"ENSG00000163584" ,"ENSG00000175467" ,"ENSG00000170779" ,"ENSG00000128228" ,"ENSG00000108344" ,"ENSG00000129474" ,"ENSG00000106615" ,"ENSG00000133030" ,"ENSG00000133265" ,"ENSG00000145494" ,"ENSG00000099956" ,"ENSG00000149925" ,"ENSG00000103196" ,"ENSG00000132475" ,"ENSG00000082153" ,"ENSG00000109685" ,"ENSG00000145741" ,"ENSG00000111642" ,"ENSG00000117632" ,"ENSG00000091527" ,"ENSG00000162302" ,"ENSG00000168883" ,"ENSG00000137166" ,"ENSG00000113140" ,"ENSG00000179950" ,"ENSG00000198053" ,"ENSG00000116260" ,"ENSG00000130164" ,"ENSG00000147140" ,"ENSG00000001617" ,"ENSG00000089685" ,"ENSG00000074410" ,"ENSG00000138180" ,"ENSG00000184007" ,"ENSG00000077348" ,"ENSG00000164099" ,"ENSG00000189091" ,"ENSG00000123080" ,"ENSG00000087077" ,"ENSG00000134291" ,"ENSG00000205250" ,"ENSG00000084623" ,"ENSG00000170955" ,"ENSG00000148773" ,"ENSG00000167468" ,"ENSG00000198055" ,"ENSG00000111674" ,"ENSG00000165943" ,"ENSG00000173530" ,"ENSG00000125901" ,"ENSG00000168994" ,"ENSG00000126461" ,"ENSG00000131504" ,"ENSG00000137804" ,"ENSG00000115306" ,"ENSG00000165119" ,"ENSG00000138326" ,"ENSG00000084774" ,"ENSG00000117984" ,"ENSG00000071553" ,"ENSG00000047849" ,"ENSG00000174136" ,"ENSG00000105401" ,"ENSG00000111057" ,"ENSG00000237493" ,"ENSG00000133112" ,"ENSG00000139433" ,"ENSG00000114631" ,"ENSG00000124541" ,"ENSG00000156976" ,"ENSG00000068394" ,"ENSG00000180900" ,"ENSG00000174684" ,"ENSG00000234664" ,"ENSG00000177731" ,"ENSG00000162520" ,"ENSG00000094975" ,"ENSG00000204381" ,"ENSG00000239306" ,"ENSG00000131711" ,"ENSG00000173113" ,"ENSG00000147872" ,"ENSG00000182095" ,"ENSG00000264364" ,"ENSG00000160957" ,"ENSG00000163814" ,"ENSG00000085840" ,"ENSG00000134853" ,"ENSG00000172667" ,"ENSG00000169564" ,"ENSG00000197930" ,"ENSG00000103888" ,"ENSG00000137054" ,"ENSG00000123485" ,"ENSG00000071655" ,"ENSG00000156599" ,"ENSG00000196498" ,"ENSG00000141682" ,"ENSG00000104408" ,"ENSG00000179348" ,"ENSG00000134668" ,"ENSG00000162576" ,"ENSG00000124733" ,"ENSG00000164087" ,"ENSG00000171148" ,"ENSG00000204174" ,"ENSG00000161642" ,"ENSG00000130517" ,"ENSG00000124731" ,"ENSG00000116688" ,"ENSG00000163931" ,"ENSG00000136068" ,"ENSG00000237506" ,"ENSG00000114019" ,"ENSG00000063046" ,"ENSG00000141556" ,"ENSG00000188760" ,"ENSG00000160256" ,"ENSG00000106636" ,"ENSG00000141560" ,"ENSG00000142002" ,"ENSG00000147955" ,"ENSG00000135506" ,"ENSG00000100823" ,"ENSG00000173706" ,"ENSG00000116288" ,"ENSG00000138685" ,"ENSG00000071564" ,"ENSG00000178952" ,"ENSG00000107854" ,"ENSG00000159231" ,"ENSG00000163110" ,"ENSG00000198763" ,"ENSG00000230202" ,"ENSG00000072958" ,"ENSG00000197905" ,"ENSG00000100138" ,"ENSG00000179958" ,"ENSG00000127955" ,"ENSG00000092470" ,"ENSG00000169180" ,"ENSG00000117411" ,"ENSG00000111665" ,"ENSG00000172009" ,"ENSG00000100216" ,"ENSG00000154473" ,"ENSG00000145901" ,"ENSG00000134107" ,"ENSG00000105835" ,"ENSG00000204176" ,"ENSG00000163703" ,"ENSG00000164615" ,"ENSG00000146701" ,"ENSG00000172594" ,"ENSG00000072163" ,"ENSG00000112893" ,"ENSG00000147604" ,"ENSG00000161847" ,"ENSG00000104824" ,"ENSG00000205581" ,"ENSG00000170004" ,"ENSG00000169714" ,"ENSG00000142156" ,"ENSG00000188643" ,"ENSG00000030582" ,"ENSG00000130811" ,"ENSG00000184967" ,"ENSG00000143382" ,"ENSG00000179750" ,"ENSG00000125356" ,"ENSG00000085733" ,"ENSG00000119335" ,"ENSG00000100034" ,"ENSG00000124225" ,"ENSG00000161888" ,"ENSG00000252645" ,"ENSG00000065911" ,"ENSG00000259918" ,"ENSG00000148908" ,"ENSG00000100528" ,"ENSG00000118526" ,"ENSG00000120937" ,"ENSG00000184564" ,"ENSG00000115241" ,"ENSG00000165271" ,"ENSG00000182871" ,"ENSG00000068438" ,"ENSG00000143320" ,"ENSG00000153187" ,"ENSG00000090889" ,"ENSG00000139278" ,"ENSG00000136153" ,"ENSG00000182827" ,"ENSG00000162613" ,"ENSG00000113657" ,"ENSG00000163923" ,"ENSG00000172336" ,"ENSG00000167085" ,"ENSG00000088986" ,"ENSG00000166913" ,"ENSG00000160075" ,"ENSG00000143612" ,"ENSG00000175792" ,"ENSG00000100285" ,"ENSG00000075223" ,"ENSG00000182985" ,"ENSG00000130429" ,"ENSG00000213846" ,"ENSG00000163661" ,"ENSG00000070087" ,"ENSG00000122378" ,"ENSG00000112578" ,"ENSG00000170855" ,"ENSG00000125877" ,"ENSG00000198682" ,"ENSG00000117298" ,"ENSG00000131747" ,"ENSG00000170017" ,"ENSG00000175602" ,"ENSG00000102898" ,"ENSG00000143314" ,"ENSG00000092969" ,"ENSG00000112139" ,"ENSG00000137936" ,"ENSG00000196549" ,"ENSG00000143256" ,"ENSG00000221829" ,"ENSG00000099331" ,"ENSG00000180964" ,"ENSG00000183431" ,"ENSG00000108312" ,"ENSG00000228253" ,"ENSG00000161016" ,"ENSG00000011028" ,"ENSG00000087365" ,"ENSG00000104341" ,"ENSG00000143369" ,"ENSG00000213190" ,"ENSG00000100591" ,"ENSG00000261295" ,"ENSG00000251562" ,"ENSG00000106771" ,"ENSG00000156711" ,"ENSG00000074696" ,"ENSG00000242265" ,"ENSG00000196365" ,"ENSG00000071626" ,"ENSG00000141367" ,"ENSG00000011426" ,"ENSG00000140988" ,"ENSG00000129038" ,"ENSG00000105486" ,"ENSG00000106333" ,"ENSG00000100764" ,"ENSG00000110090" ,"ENSG00000105197" ,"ENSG00000143622" ,"ENSG00000197321" ,"ENSG00000188986" ,"ENSG00000109572" ,"ENSG00000112559" ,"ENSG00000186283" ,"ENSG00000272327" ,"ENSG00000123143" ,"ENSG00000139289" ,"ENSG00000198825" ,"ENSG00000064601" ,"ENSG00000169567" ,"ENSG00000173456" ,"ENSG00000221869" ,"ENSG00000168216" ,"ENSG00000169991" ,"ENSG00000106546" ,"ENSG00000143476" ,"ENSG00000089159" ,"ENSG00000132383" ,"ENSG00000079246" ,"ENSG00000065135" ,"ENSG00000232480" ,"ENSG00000122566" ,"ENSG00000159259" ,"ENSG00000142798" ,"ENSG00000146535" ,"ENSG00000116670" ,"ENSG00000185022" ,"ENSG00000155506" ,"ENSG00000167747" ,"ENSG00000158186" ,"ENSG00000130826" ,"ENSG00000184207" ,"ENSG00000071282" ,"ENSG00000007520" ,"ENSG00000084754" ,"ENSG00000135446" ,"ENSG00000134686" ,"ENSG00000144381" ,"ENSG00000130332" ,"ENSG00000213024" ,"ENSG00000197771" ,"ENSG00000224861" ,"ENSG00000007923" ,"ENSG00000206503" ,"ENSG00000126698" ,"ENSG00000065000" ,"ENSG00000107223" ,"ENSG00000165434" ,"ENSG00000101460" ,"ENSG00000141753" ,"ENSG00000089154" ,"ENSG00000101365" ,"ENSG00000247077" ,"ENSG00000228305" ,"ENSG00000166503" ,"ENSG00000138735" ,"ENSG00000084207" ,"ENSG00000104332" ,"ENSG00000203485" ,"ENSG00000242372" ,"ENSG00000130520" ,"ENSG00000100196" ,"ENSG00000136930" ,"ENSG00000071082" ,"ENSG00000101849" ,"ENSG00000142552" ,"ENSG00000177954" ,"ENSG00000153714" ,"ENSG00000148516" ,"ENSG00000169813" ,"ENSG00000166482" ,"ENSG00000108932" ,"ENSG00000130224" ,"ENSG00000006468" ,"ENSG00000232573" ,"ENSG00000163453" ,"ENSG00000115226" ,"ENSG00000091483" ,"ENSG00000151892" ,"ENSG00000145050" ,"ENSG00000166510" ,"ENSG00000063978" ,"ENSG00000142634" ,"ENSG00000138448" ,"ENSG00000135074" ,"ENSG00000119138" ,"ENSG00000130956" ,"ENSG00000146918" ,"ENSG00000116044" ,"ENSG00000122687" ,"ENSG00000149571" ,"ENSG00000091542" ,"ENSG00000105438" ,"ENSG00000244586" ,"ENSG00000197702" ,"ENSG00000060688" ,"ENSG00000148730" ,"ENSG00000184887" ,"ENSG00000123064" ,"ENSG00000044574" ,"ENSG00000134871" ,"ENSG00000136205" ,"ENSG00000064961" ,"ENSG00000178921" ,"ENSG00000146072" ,"ENSG00000104419" ,"ENSG00000132031" ,"ENSG00000122863" ,"ENSG00000174695" ,"ENSG00000090971" ,"ENSG00000159228" ,"ENSG00000198814" ,"ENSG00000149480" ,"ENSG00000008441" ,"ENSG00000065882" ,"ENSG00000165655" ,"ENSG00000189306" ,"ENSG00000204253" ,"ENSG00000154978" ,"ENSG00000166557" ,"ENSG00000174444" ,"ENSG00000069849" ,"ENSG00000143753" ,"ENSG00000188460" ,"ENSG00000273179" ,"ENSG00000065613" ,"ENSG00000138593" ,"ENSG00000186340" ,"ENSG00000177606" ,"ENSG00000109814" ,"ENSG00000148841" ,"ENSG00000141522" ,"ENSG00000004455" ,"ENSG00000022277" ,"ENSG00000253304" ,"ENSG00000170727" ,"ENSG00000259781" ,"ENSG00000164182" ,"ENSG00000136542" ,"ENSG00000116157" ,"ENSG00000167693" ,"ENSG00000132361" ,"ENSG00000102900" ,"ENSG00000160087" ,"ENSG00000112245" ,"ENSG00000143515" ,"ENSG00000138829" ,"ENSG00000124098" ,"ENSG00000106004" ,"ENSG00000113645" ,"ENSG00000160200" ,"ENSG00000185633" ,"ENSG00000132824" ,"ENSG00000071859" ,"ENSG00000140105" ,"ENSG00000167325" ,"ENSG00000065268" ,"ENSG00000257219" ,"ENSG00000132002" ,"ENSG00000077942" ,"ENSG00000089693" ,"ENSG00000167657" ,"ENSG00000107833" ,"ENSG00000031081" ,"ENSG00000125538" ,"ENSG00000029153" ,"ENSG00000174080" ,"ENSG00000144583" ,"ENSG00000125821" ,"ENSG00000007255" ,"ENSG00000120802" ,"ENSG00000143436" ,"ENSG00000130766" ,"ENSG00000100726" ,"ENSG00000197226" ,"ENSG00000147274" ,"ENSG00000198952" ,"ENSG00000099995" ,"ENSG00000134308" ,"ENSG00000137713" ,"ENSG00000150281" ,"ENSG00000008256" ,"ENSG00000116285" ,"ENSG00000159377" ,"ENSG00000101811" ,"ENSG00000051180" ,"ENSG00000180875" ,"ENSG00000268942" ,"ENSG00000143067" ,"ENSG00000154146" ,"ENSG00000127616" ,"ENSG00000104812" ,"ENSG00000261061" ,"ENSG00000134809" ,"ENSG00000122034" ,"ENSG00000116489" ,"ENSG00000168758" ,"ENSG00000100441" ,"ENSG00000111716" ,"ENSG00000106803" ,"ENSG00000168268" ,"ENSG00000183779" ,"ENSG00000142208" ,"ENSG00000099875" ,"ENSG00000114686" ,"ENSG00000171067" ,"ENSG00000167797" ,"ENSG00000171241" ,"ENSG00000068912" ,"ENSG00000187778" ,"ENSG00000077238" ,"ENSG00000063245" ,"ENSG00000130202" ,"ENSG00000114520" ,"ENSG00000168610" ,"ENSG00000119403" ,"ENSG00000106290" ,"ENSG00000003056" ,"ENSG00000179271" ,"ENSG00000123159" ,"ENSG00000163191" ,"ENSG00000196363" ,"ENSG00000163083" ,"ENSG00000214182" ,"ENSG00000138678" ,"ENSG00000189184" ,"ENSG00000166813" ,"ENSG00000124299" ,"ENSG00000148335" ,"ENSG00000105447" ,"ENSG00000176845" ,"ENSG00000134779" ,"ENSG00000160113" ,"ENSG00000111885" ,"ENSG00000114023" ,"ENSG00000133026" ,"ENSG00000172534" ,"ENSG00000156508" ,"ENSG00000070614" ,"ENSG00000131791" ,"ENSG00000178464" ,"ENSG00000163359" ,"ENSG00000189343" ,"ENSG00000164574" ,"ENSG00000116741" ,"ENSG00000085999" ,"ENSG00000168300" ,"ENSG00000182446" ,"ENSG00000179051" ,"ENSG00000175756" ,"ENSG00000139211" ,"ENSG00000119899" ,"ENSG00000227081" ,"ENSG00000067167" ,"ENSG00000203813" ,"ENSG00000119541" ,"ENSG00000008952" ,"ENSG00000132003" ,"ENSG00000119280" ,"ENSG00000142541" ,"ENSG00000006744" ,"ENSG00000162063" ,"ENSG00000113719" ,"ENSG00000197256" ,"ENSG00000226415" ,"ENSG00000099817" ,"ENSG00000138834" ,"ENSG00000182512" ,"ENSG00000162244" ,"ENSG00000111275" ,"ENSG00000071967" ,"ENSG00000100979" ,"ENSG00000262814" ,"ENSG00000198873" ,"ENSG00000135829" ,"ENSG00000113013" ,"ENSG00000164587" ,"ENSG00000213412" ,"ENSG00000111786" ,"ENSG00000153113" ,"ENSG00000089737" ,"ENSG00000133639" ,"ENSG00000185651" ,"ENSG00000101224" ,"ENSG00000198860" ,"ENSG00000196396" ,"ENSG00000101439" ,"ENSG00000184838" ,"ENSG00000138623" ,"ENSG00000167881" ,"ENSG00000095319" ,"ENSG00000138629" ,"ENSG00000164040" ,"ENSG00000168005" ,"ENSG00000107165" ,"ENSG00000170373" ,"ENSG00000163382" ,"ENSG00000023191" ,"ENSG00000111737" ,"ENSG00000105722" ,"ENSG00000108270" ,"ENSG00000125503" ,"ENSG00000151176" ,"ENSG00000143569" ,"ENSG00000049449" ,"ENSG00000108774" ,"ENSG00000154319" ,"ENSG00000090621" ,"ENSG00000074416" ,"ENSG00000176393" ,"ENSG00000135451" ,"ENSG00000113578" ,"ENSG00000168874" ,"ENSG00000000003" ,"ENSG00000108960" ,"ENSG00000141985" ,"ENSG00000134824" ,"ENSG00000038427" ,"ENSG00000137193" ,"ENSG00000162627" ,"ENSG00000197956" ,"ENSG00000166311" ,"ENSG00000197632" ,"ENSG00000076248" ,"ENSG00000219201" ,"ENSG00000180596" ,"ENSG00000058262"]
    # gene_list =gene_list[:100]
    network_file_name = "PCNet.sif"
    truncate_network_by_gene_list(gene_list, network_file_name, radius=5)







