#!/usr/bin/ruby

#require "rubygems"
#require "intermine/service"
#service = Service.new("https://apps.araport.org:443/thalemine")

# Get a new query from the service you will be querying:


require "rubygems"
require "intermine/service"
service = Service.new("https://apps.araport.org:443/thalemine")
service.new_query("Protein").
    select(["name", "primaryAccession", "genes.primaryIdentifier", "genes.symbol", "genes.briefDescription", "genes.isObsolete"]).
    order_by("name", "ASC").
#    limit(10).
    each_row do |r| 
      if r["genes.primaryIdentifier"] =~ /^AT\dG/ 
        puts r["genes.primaryIdentifier"]
      end
    end


#service.new_query("Protein").
#    select(["name", "primaryAccession", "genes.primaryIdentifier", "genes.symbol", "genes.briefDescription", "genes.isObsolete"]).
    # You can edit the constraint values below
#    where("genes.isObsolete" => {"==" => false}).
#    order_by("primaryIdentifier", "ASC").
#    each_row do |r| 
#      if r["genes.primaryIdentifier"] =~ /^AT\dG/ 
#        puts r["genes.primaryIdentifier"]
#      end
#    end