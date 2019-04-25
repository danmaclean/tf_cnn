directory "lib"

file "lib/affy_ATH1_array_elements-2010-12-20.txt" => "lib" do
  sh %{wget -P lib/ https://www.arabidopsis.org/download_files/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt }
end

file "lib/AtRegNet" => "lib" do
  sh %{wget -P lib/ http://agris-knowledgebase.org/Downloads/AtRegNet.zip }
end

file "lib/Ath_TF_list" => "lib" do
  sh %{wget -P lib/ http://planttfdb.cbi.pku.edu.cn/download/TF_list/Ath_TF_list.gz }
end

desc "gets agis for all proteins in Araport"
file "lib/all_agis.txt" => "lib" do
  sh %{ruby scripts/get_gene_names.rb > lib/all_agis.txt }
end

desc "runs predictions for all TF/Target pairs"
file "data/predictions.csv" => "data" do
  sh %{r i in `seq 1 1717`; do Rscript scripts/run_single_tf_predict.R $i done}
  sh %{cat tmp/*.csv data/predictions.csv}
end

task :prepare_data do
  Rake::Task["lib/affy_ATH1_array_elements-2010-12-20.txt"].invoke
  Rake::Task["lib/AtRegNet"].invoke
  Rake::Task["lib/all_agis.txt"].invoke
  Rake::Task["data/arab_TF.RDS"].invoke
  Rake::Task["data/scaled_arab_TF.RDS"].invoke
end

task "data/arab_TF.RDS" => ["data", "lib/AtRegNet", "lib/affy_ATH1_array_elements-2010-12-20.txt"] do
  sh %{Rscript scripts/make_training_data.R}
end

task "data/scaled_arab_TF.RDS" => ["data", "lib/AtRegNet", "lib/affy_ATH1_array_elements-2010-12-20.txt"] do
  sh %{Rscript scripts/make_training_data.R scale}
end


task "data/scaled_arab_TF_two_channel.RDS" => ["data", "lib/AtRegNet", "lib/affy_ATH1_array_elements-2010-12-20.txt"] do
  sh %{Rscript scripts/make_training_data_for_convnet.R scale}
end

task "data/scaled_arab_TF_one_channel.RDS" => ["data", "lib/AtRegNet", "lib/affy_ATH1_array_elements-2010-12-20.txt"] do
  sh %{Rscript scripts/make_training_data_for_convnet_one_channel.R scale}
end