directory "lib"

file "lib/affy_ATH1_array_elements-2010-12-20.txt" => "lib" do
  sh %{wget -P lib/ https://www.arabidopsis.org/download_files/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt }
end

file "lib/AtRegNet" => "lib" do
  sh %{wget -P lib/ http://agris-knowledgebase.org/Downloads/AtRegNet.zip }
end

task :prepare_data do
  Rake::Task["lib/affy_ATH1_array_elements-2010-12-20.txt"].invoke
  Rake::Task["lib/AtRegNet"].invoke
  Rake::Task["data/arab_TF.RDS"]
  Rake::Task["data/scaled_arab_TF.RDS"]
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