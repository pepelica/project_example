# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:30:52 2020

@author: Aleksandra
"""


import os
from argparse import ArgumentParser
import shutil
import time
import functions as fu
import configparser


quantile = 4
 
def main(quantile):
    
    parser = ArgumentParser(description = "Find quantiles for proteins expression in gbk files.")
    
    parser.add_argument('-i', '--indir', type = str,
                        help='A directory containing folders with gbk files.')
    parser.add_argument('-o', '--outdir', type = str,
                        help="A directory for result files.")
    parser.add_argument('-q', '--quantiles', type = int, default = 4,
                        help="Number of quartiles (default: %(default)s).")


    args = parser.parse_args()
    print(args)
    
    #read ini file
    config = configparser.ConfigParser() 
    ini_path = os.path.dirname(os.path.realpath(__file__))
    config.read(ini_path + '\path.ini')
    eloe_path = config['PATH']['eloe_path']
    input_ini_path = config['PATH']['input_path']
    output_ini_path = config['PATH']['output_path']
    
    
    input_path = args.indir
    output_path = args.outdir
    quantile = args.quantiles
    
    #если пути не введены в командную строку, берутся пути из ini файла
    if not input_path:
        input_path = input_ini_path
    if not output_path:
        output_path = output_ini_path

    try:                                    
        samples = os.listdir(input_path)   #проверка, что входной путь корректен
    except(FileNotFoundError):
        print('Please provide a command line parameter using --indir or -i or add input path to the path.ini file')
        return
    
    #удаляем все файлы в рабочей директории
    try:                                    #проверка, что выходной путь корректен
        contents = os.listdir(output_path)
    except(FileNotFoundError):
        print('Please provide a command line parameter using --outdir or -o or add output path to the path.ini file')
        return
    
    if not os.path.exists(eloe_path):
        print('Please add correct eloe path to the path.ini file')
        return
        
    if len(contents)!= 0:
        for file in contents:
            shutil.rmtree(output_path + '\\'+file)
            
    
    #Создаем общую папку для результата подсчета кванитлей
    quantile_output_path = output_path+'\\quantile\\'
     
    os.mkdir(quantile_output_path)
    


    for sample in samples:
        
        path = input_path+ '\\' + sample
        org_temp_res_path = output_path + '\\temp\\'+ sample
        
        
        gbk_files = org_temp_res_path + '\\gbk\\'+sample
        try:
            os.makedirs(gbk_files)
        except:
            shutil.rmtree(output_path + '\\temp\\'+sample)
            os.makedirs(gbk_files)  
            
        try:

            #читаем файл, если некорректен - преобразуем, сохраняем в новый файл
            #возвращает список геномных карточек
            records = fu.read_input_file(path, gbk_files+'\\normal.gbk') #read and modificate file for parsing
      

        except ValueError:
            print(sample + ' is empty')
                  
        except:       
            print(sample +' have incorrect file format')
        
        else:
            #cклеиваем геномные карточки в одну, записываем в файл
            fu.merge(records, gbk_files+'\\merged.gbk') #merge records to one record 

                
            #если исправлялся некорректный файл - удаляем полученный файл
            #(чтобы в папке, подаваемой елое был только один файл)
            try:
                os.remove(gbk_files+'\\normal.gbk')
            except:
                print('file was not modified')
             
            #создаем папку с результатами елое этого организма
            eloe_res_path = org_temp_res_path+'\\eloe\\'
            os.makedirs(eloe_res_path)
             
            #запускаем елое, записываем результаты для организма
            eloe_res = fu.run_EloE(org_temp_res_path + '\\gbk', eloe_res_path, eloe_path)
            
            #папка, куда записываются файлы с квантилями этого организма 
            sample_quantile = quantile_output_path + '\\'+ sample+'\\'
            os.mkdir(sample_quantile) 
            if eloe_res:
                fu.find_quantile(eloe_res_path+sample.replace(' ', '_'), sample_quantile + sample, quantile)
                print('Found quantile for sample '+sample)
            else: 
                print ("EloE can't finish the process for "+sample)
            
        
       

    
    #удаляет рабочую папку с gbk файлами и результатами eloe
    shutil.rmtree(output_path + '\\temp')
        


main(quantile)