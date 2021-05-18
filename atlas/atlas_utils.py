import os
import pandas as pd
import multiprocessing
import subprocess
from pathlib import Path


class Utils():

    def __init__(self, na, nr, rmax, desc, output_path, bof):

        self.bco = bof
        self.pc = "../../../../bin/bin/phase-congruency"
        self.land = "../../../../bin/bin/landmark-detector"
        self.prob_atlas = "../../../../bin/bin/probabilistic-landmark-atlas"
        self.output = output_path
        self.na, self.nr, self.rmax = [na, nr, rmax]
        self.desc = desc

# -------------------------------------------------------------------------------------------------
# Detection of PC landmarks for each input image
# -------------------------------------------------------------------------------------------------

    def multip(self, input1, input2, path_out_, func, n_proc):
        params_ = zip(input1, input2, path_out_)
        p = multiprocessing.Pool(n_proc)
        p.map(func, params_)
        p.close()

    def get_directories_mask(self, root_folder, root_mask, new_folder):
        list_in, list_mask, list_out = [], [], []
        for filename in os.listdir(root_folder):
            if "_mask.nii.gz" in filename:
                filename = filename.split("_mask.nii.gz")[0]
                list_mask.append(root_mask + "/" + filename)
                list_in.append(root_folder + "/" + filename)
                list_out.append(new_folder + "/" + filename)
        return list_in, list_mask, list_out

    def phase_congruency(self, params_):
        inputfile, _, outputfile = params_
        if not Path(outputfile+"_eigenvalues_0.nii.gz").exists():
            print("time", self.pc, "-i", inputfile, "-bof", self.bco,
                  "-o", outputfile)
            subprocess.call(["time", self.pc, "-i", inputfile,
                             "-bof", self.bco, "-o", outputfile])

    def landmarks(self, params_):
        inputfile, input_mask, inputpc = params_
        output_name = "/" + inputpc.split("/")[-1]
        if not Path(self.output + "/" + output_name +
                    "_landmarks.txt").exists():
            print("time", self.land, "-r", inputfile, "-io", inputpc,
                  self.desc, "-st", "0.33", "-nr", self.nr,
                  "-na", self.na, "-gl", self.rmax, "-getp",
                  "-mask", input_mask, "-o", self.output + output_name)
            subprocess.call(["time", self.land, "-r", inputfile, "-io",
                             inputpc, self.desc, "-st", "0.33",
                             "-nr", self.nr, "-na", self.na, "-gl",
                             self.rmax, "-getp", "-mask", input_mask,
                             "-o", self.output + output_name])

    def create_folders(self, new_folder):
        if not os.path.exists(new_folder):
            os.makedirs(new_folder, exist_ok=True)

# -------------------------------------------------------------------------------------------------
# Construction of the average and probabilistic atlases
# -------------------------------------------------------------------------------------------------

    def get_directories_landmarks(self, root_folder):
        list_in = []
        for name in os.listdir(root_folder):
            if "_landmarks.txt" in name:
                list_in.append(name)
        return list_in

    def get_input_paths_atlas(self, output_age_path, hospital_names,
                              auxdata_path, atlas_samples_prefix):
        inputs = self.get_directories_landmarks("output")

        self.create_folders(output_age_path)
        m = 0

        reader = pd.read_csv(auxdata_path+"/IXI_20-80.txt",
                             header=None, dtype=str)
        for hospital in hospital_names:
            output_age_hosp_path = output_age_path+"/"+hospital
            self.create_folders(output_age_hosp_path)
            for i in reader.index:
                subject_id = reader.iloc[i].values[0]
                img_name = "IXI"+str(subject_id)+"-"+hospital

                res = [v for v, in_ in enumerate(inputs) if img_name in in_]
                if len(res) > 0:
                    os.system("cp output/" + inputs[res[0]] + " " +
                              output_age_path + "/" + atlas_samples_prefix +
                              "_" + str(m) + ".txt")
                    m += 1
        return m

    def probabilistic_atlas_landmarks(self, output_age_path, m,
                                      atlas_samples_prefix, template_image):
        print("time", self.prob_atlas, "-i",
              output_age_path+"/"+atlas_samples_prefix,
              "-n", str(m), "-r", template_image, "-kstd", "1.0",
              "-rstd", "1.0", "-ek", "50", "-ck", "5", "-dt",
              "0.1", "-aw", "3", "-o", "probabilistic_atlas")
        subprocess.call(["time", self.prob_atlas, "-i",
                         output_age_path+"/"+atlas_samples_prefix,
                         "-n", str(m), "-r", template_image, "-kstd",
                         "1.0", "-rstd", "1.0", "-ek", "50", "-ck", "5",
                         "-dt", "0.1", "-aw", "3", "-o",
                         "probabilistic_atlas"])
        os.system("mv probabilistic_atlas* "+output_age_path)
