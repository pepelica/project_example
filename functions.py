from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, Alphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd
from subprocess import Popen, PIPE
import os





def merge(records, path):
    """Merge multiple SeqRecords into one, using a defined spacer

    :param records: Iterable containing SeqRecords to be merged
    :param length: Length of the spacer in kbp
    :param spacer: Kind of spacer to use ('n' for UnknownSeq spacer, 'stop' for all-frame stop codon spacer)

    :return: A single SeqRecord that is the product of the merge.
    """
    length=20
    spacer='n'
    ALL_FRAME_STOP_MOTIF = 'TAGCTAACTGACCGTCAGTTAGCTA'

    if not len(records):
        raise ValueError("No records given")

    if spacer == 'stop':
        spacer_seq = Seq(ALL_FRAME_STOP_MOTIF * 40 *
                         length, Alphabet.generic_dna)
    else:
        spacer_seq = UnknownSeq(
            length * 1000, alphabet=Alphabet.generic_dna, character='N')

    new_rec = records[0]

    if len(records) == 1:
        with open(path, 'w') as merged_file:
           SeqIO.write(new_rec, merged_file, 'genbank')
    else:  
        rec_id = new_rec.id
        rec_name = new_rec.name
        rec_desc = new_rec.description
        date = new_rec.annotations.get('date', '')
        source = new_rec.annotations.get("source", '')
        organism = new_rec.annotations.get('organism', '')
        taxonomy = new_rec.annotations.get('taxonomy', [])
        data_file_division = new_rec.annotations.get('data_file_division', 'UNK')
        topology = new_rec.annotations.get('topology', 'linear')
    
        for i, rec in enumerate(records[1:]):
            spacer_id = 'spacer_{}'.format(i + 1)
    
            spacer_feature = SeqFeature(FeatureLocation(0, length * 1000, 0),
                                        type='misc_feature', id=spacer_id,
                                        qualifiers={'note': [spacer_id]})
    
            spacer_rec = SeqRecord(spacer_seq, id=spacer_id, name=spacer_id,
                                   description=spacer_id, features=[spacer_feature])
    
            new_rec = new_rec + spacer_rec + rec
    
        new_rec.id = rec_id
        new_rec.name = rec_name
        new_rec.description = rec_desc
        new_rec.annotations["date"] = date
        new_rec.annotations["source"] = source
        new_rec.annotations["organism"] = organism
        new_rec.annotations["taxonomy"] = taxonomy
        new_rec.annotations["data_file_division"] = data_file_division
        new_rec.annotations["topology"] = topology
        
        with open(path, 'w') as merged_file:
            SeqIO.write(new_rec, merged_file, 'genbank')



#сюда попадаем, если файл не парсится; создается читаемая версия файла
def modify_input_file(handle):
    def removeTrailingZeros(line):
        res=''
        ss=line.split('.')
        if len(ss)>1:
            tmp = ss[1].strip()
            tmprev = tmp[::-1]#reverse
            #читаем дробную часть в обратном порядке
            flag=True
            resrev=''
            for c in tmprev: 
                if (c=='0') & flag:
                    continue
                else:
                    flag=False
                    resrev = resrev + str(c)
            #если не было ничего кроме нулей после точке, как в примере
            if resrev=='':
                res='0'
            else:
                res=resrev[::-1]
        return ss[0]+'.'+res


    new_handle = []
    new_line = ''

    for line in handle:
        new_line = line
        if 'LOCUS       ' in line:
            if line.count('-')<2:
                line = line.replace('\n', '       ')
            len_start = line.find('length_')+7
            len_end = line.find('_cov')
            lenght = line[len_start:len_end]
            indx = line.rfind(lenght)
            head2 = ' '+ line[indx:] # вторая половина заголовка, дальше понадобится
            line = line[:indx]+head2


            if len(line.split()[1]) > 16:

                line = line.replace("NODE", "N")
                line = line.replace("length", "l")
                line = line.replace("cov", "c")

                head1 = 'LOCUS       ' + removeTrailingZeros(line.split()[1])
                line = head1 + head2
                #если и это не помогло, то убираем инфу о покрытии



                if len(line.split()[1]) > 16:
                    head1 = line.split("_c")[0]
                    line = head1 + head2
                    #если и это не помогло, то убираем инфу о длине

                    if len(line.split()[1]) > 16:
                        head = line.split("_l")[0]
                        line = head1 + head2
                        #если и это не помогло, то убираем разделители _ в названии


                        if len(line.split()[1]) > 16:
                            line = line.replace('_', '')

                            if len(line.split()[1]) > 16:
                                line = 'LOCUS       '+ line.split()[1][:16] + head2

            new_line = line



        new_handle.append(new_line)

    return ''.join(new_handle)

        
                
#считывает входной файл, если не парсится - отправляет на корректировку, пересохраняет и парсит, возвращает список геномных карточек образца                
def read_input_file(path, new_handle_file):

    contents = os.listdir(path)

    records = []
    if len(contents) == 1:

        try:
            with open(path+'//'+ contents[0], "r") as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    records.append(record)
        except:
            with open(path+'//'+ contents[0], "r") as handle:
                new_handle = modify_input_file(handle)
                with open(new_handle_file, 'w') as myfile:
                    myfile.write(new_handle)


            for record in SeqIO.parse(new_handle_file, "genbank"):
                records.append(record)
    elif not len(contents):
        raise ValueError
                    
    else:
        for genome in contents:
            with open(path + '//'+genome, 'rU') as handle:
                record = SeqIO.parse(handle, 'genbank')
                records.append(record)
    return records

#запускает программу eloe
def run_EloE(EloE_input_path, EloE_output_path, eloe_path):

    proc = Popen(
        ['java', '-jar', eloe_path,  
         EloE_input_path, EloE_output_path],
        shell=True,
        stdout=PIPE, stderr=PIPE)
    proc.wait()    # дождаться выполнения

    try:
        proc.communicate(timeout=15)       
        return True
    except Popen.TimeoutExpired:
        proc.kill()
        return False
    except:

        return False
        

#находит файл eei и разбивает его на заданное количество частей
def find_quantile(EloE_output_path, output, quantile):
    EloEresults = os.listdir(EloE_output_path)
    for EloE_file in EloEresults:
        if '_eei' in EloE_file:
            df = pd.read_csv(EloE_output_path+'\\'+EloE_file, delimiter='\t', 
                                   usecols = ['locus_tag', 'protein_id', 'gene_id',
                                              'gene_name', 'COG', 'EEI'])
            start = 0
            end = 0
            for i in range(quantile):
                start = end
                end = len(df)-round(len(df)/4*(quantile - 1-i))
#                df1 = df.iloc[round(i*len(df)/4):len(df)-round(len(df)/4*(quantile - 1-i)), :]
                df1 = df.iloc[start:end, :]

                df1.to_csv(f'{output}_quantile{i+1}.txt', sep = '\t', na_rep = '-')
            break




