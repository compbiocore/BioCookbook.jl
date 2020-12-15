
function find_pos_val(val, input_string)
    t= split(split(input_string)[8],";")
    ind = [occursin(val,i) for i in t]
    y= replace(t[ind][1],val => "")
    return y
end

"""
    get_sv_attrs(in_str::String)::Dict

Parse a line from the vcf file and return a dictionary with the following attributes as keys
    - start: SV start position
    - end: SV end position
    - chr: SV Chromosome
    - svlen: Length tof the SV
    - svtype: The SV type

"""
function get_sv_attrs(in_str)
    dict = Dict()
    dict["chr"] = split(in_str)[1]
    dict["start"] = split(in_str)[2]
    dict["end"] = find_pos_val("END=",in_str)
    dict["svlen"] = find_pos_val("SVLEN=",in_str)
    dict["svtype"] = find_pos_val("SVTYPE=",in_str)
    dict["seq"] = split(in_str)[5]
    return dict
end


function annotate_sv(sv,x, y)
    if sv =="Unique" || sv == "unclassified"
        return x*"\t"*sv*"\tNA\tNA"

    else
        return x*"\t"*sv*"\t"*replace(y,"\t" => "|") 
    end
end


"""
    compare_locations(x::String, y:String)::String

Compare two lines each from a single vcf file and return what the status of the location is. The function checks for the following 
  - Unique: start and end not equal
  - Merged: Start and end equal and of same SV type
  - singleton: Either one of start or end matches and of the same SV type
  - unclassified: All other cases
"""
function compare_locations(x,y)
    x_attrs = get_sv_attrs(x)
    y_attrs = get_sv_attrs(y)
    if x_attrs["chr"] != y_attrs["chr"]
        return "Unique"
    end 
    #println(x)
    consensus = x_attrs["seq"] == y_attrs["seq"]
    if x_attrs["start"] != y_attrs["start"] &&  x_attrs["end"] != y_attrs["end"] 
        return "Unique"

    elseif x_attrs["start"]== y_attrs["start"] && x_attrs["end"] == y_attrs["end"] && x_attrs["svtype"] == y_attrs["svtype"]
        return "Both"*"\t"*string(consensus)

    elseif  x_attrs["start"]== y_attrs["start"] && x_attrs["end"] != y_attrs["end"] && x_attrs["svtype"] == y_attrs["svtype"]
        return "Left"*"\t"*string(consensus)

    elseif  x_attrs["start"]!= y_attrs["start"] && x_attrs["end"] == y_attrs["end"] && x_attrs["svtype"] == y_attrs["svtype"]
        return "Right"*"\t"*string(consensus)

    else
        return "unclassified"
    end
end

"""
    check_pos(x:String, y::String)::Bool

Compare two input SV lines to make sure that the position of the first line is lower than or equal to the other
"""
function check_pos(x,y)
    chr_start_x = split(x)[1]
    pos_x = split(x)[2]
    chr_start_y = split(y[1])[1]
    pos_y = split(y[1])[2]
    if parse(Int,pos_x) <= parse(Int,pos_y) && chr_start_x == chr_start_y
        return true
    elseif parse(Int,pos_x ) >  parse(Int,pos_y) && chr_start_x != chr_start_y
        return true
    else
        return false
    end
end


function annotate_pindel(v1,v2)
    check_line = iterate(v2)
    v3=[]
    i = 1
    j = 1
    for line1 in v1
        if check_pos(line1,check_line) 
            sv = compare_locations(line1,check_line[1])
        elseif !check_pos(line1,check_line) && iterate(v2, check_line[2]) != nothing 
            while !check_pos(line1,check_line)
                check_line = iterate(v2,check_line[2])
            end
            sv= compare_locations(line1,check_line[1])
        else
            sv= compare_locations(line1,check_line[1])
        end
        val = annotate_sv(sv,line1,check_line[1])
        #println(val)
        push!(v3,val)
        #println(string(i)*" " *string(j))
        i +=1
    end
    return v3
end

function print_res(filePath,vec)
    fheader = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT\tAD_REF\tAD_ALT"
    fheader = fheader*"\tStatus_of_merge\tSeq_overlap_status\tResults_of_merge"
    fheader = fheader*"\tSV_Type\tSV_length\tHompolymer_length\tfrac_ALT"
    f = open(filePath,"w")
    println(f,fheader)
    for line in vec
        tmpline = split(line)
        tmpline[10] =replace(tmpline[10],":" => "\t")
        tmpline[10] =replace(tmpline[10],"," => "\t")
        AD_ref = parse(Float64,split(tmpline[10])[2])
        AD_alt = parse(Float64,split(tmpline[10])[3])
        
        tmpline_to_print = tmpline
        tmpline_to_print[4] = first(tmpline_to_print[4],10)
        tmpline_to_print[5] = first(tmpline_to_print[5],10)
        tmpline_to_print[13] = first(tmpline_to_print[13],10)

        push!(tmpline_to_print, find_pos_val("SVTYPE=",line))
        push!(tmpline_to_print, find_pos_val("SVLEN=",line))
        push!(tmpline_to_print, find_pos_val("HOMLEN=",line))
        push!(tmpline_to_print, string(AD_alt/(AD_alt+AD_ref)))
        println(f,join(tmpline_to_print,"\t"))
    end
    close(f)
end

function run_pindel_compare(samp1, samp2; datdir="/Users/aragaven/scratch/test_julia/test_pindel/pindel_vcf", 
                            outdir= "/Users/aragaven/Documents/Research/CBC/ene/")
    #datdir = "/Users/aragaven/scratch/test_julia/test_pindel/pindel_vcf";
    #samp1 = "DSY303_OP";
    #samp2 = "DSY303_WH";


    vec1 = readlines(joinpath(datdir,samp1*".vcf"));
    vec2 = readlines(joinpath(datdir,samp2*".vcf"));

    vec1_1 = filter(x -> !startswith(x, "#"), vec1);
    vec2_2 = filter(x -> !startswith(x, "#"), vec2);


    vec4 = annotate_pindel(vec1_1,vec2_2);
    vec5 = annotate_pindel(vec2_2,vec1_1);

    filep1 = outdir*"/"*samp1*".results.tsv"
    print_res(filep1,vec4)
    filep2 = outdir*"/"*samp2*".results.tsv"
    print_res(filep2,vec5)
end

#run_pindel_compare("DSY303_OP","DSY303_WH")
run_pindel_compare("RBY1118_OP","RBY1118_WH")
