function comp_loc(x,y, dat_col)
    chr_start = split(x)[1]
    chr_start_2 = split(y)[1]
    if chr_start != chr_start_2
        return "Unique"
    end  
    #println(x)
    pos_start = split(x)[2]
    #pos_end = replace(split(split(x)[8],";")[1],"END=" => "")
    pos_end = find_pos_val("END=",x)
    #svlen = replace(split(split(x)[8],";")[4],"SVLEN=" => "")
    svlen = find_pos_val("SVLEN=",x)
    svtype = find_pos_val("SVTYPE=",x)
    if occursin("INV",svtype)
    elseif occursin("DEL",svtype)#svtype != "DEL" 
    elseif occursin("DUP:TANDEM",svtype)
    else
        pos_end = string(parse(Int, pos_end) + parse(Int,svlen))
    end
    #=if svtype != "DEL" 
        pos_end = string(parse(Int, pos_end) + parse(Int,svlen))
    end=#

    pos_start_2 = split(split(split(split(y)[dat_col],":")[11],"-")[1],"_")[end]
    pos_end_2 = split(split(split(split(y)[dat_col],":")[11],"-")[2],"_")[end]
    if pos_start != pos_start_2 &&  pos_end != pos_end_2
        return "Unique"
    elseif pos_start == pos_start_2 && pos_end == pos_end_2
        return "Merged"
    else
        #println("error")
        #println(svtype*"\t"*pos_start*"\t"*pos_end*"\t"*pos_start_2*"\t"*pos_end_2)
        #println(x)
        return "singleton"
    end
end

function find_pos_val(val, instring)
    t= split(split(instring)[8],";")
    ind = [occursin(val,i) for i in t]
    y= replace(t[ind][1],val => "")
    return y
end

"""
    check_pos(x:String, y::String)::Bool

Compare two input SV lines to make sure that the position of the first line is lower than or equal to the other


"""
function check_pos(x,y)
    chr_start = split(x)[1]
    pos_1 = split(x)[2]
    chr_start_2 = split(y[1])[1]
    pos_2 = split(y[1])[2]
    if parse(Int,pos_1 ) <= parse(Int,pos_2) && chr_start == chr_start_2
        return true
    elseif parse(Int,pos_1 ) >  parse(Int,pos_2) && chr_start != chr_start_2
        return true
    else
        return false
    end
end

function compare_break_points(x,y,comp_col)
    pos_start_X = split(x)[2]
    #pos_end = replace(split(split(x)[8],";")[1],"END=" => "")
    pos_end_x = find_pos_val("END=",x)
    svlen_x = find_pos_val("SVLEN=",x)
    svtype_x = find_pos_val("SVTYPE=",x)
    if svtype_x != "INV"
        pos_end_x = string(parse(Int, pos_end_x) + parse(Int,svlen_x))
    end

    pos_start_y = split(split(split(split(y)[comp_col],":")[11],"-")[1],"_")[end]
    pos_end_y = split(split(split(split(y)[comp_col],":")[11],"-")[2],"_")[end]
    if pos_start_x == pos_start_y && pos_end_x == pos_end_y
        bk_pt_overlap ="Both"
    elseif pos_start_x == pos_start_y && pos_end_x != pos_end_y
        bk_pt_overlap = "Left"
    elseif pos_end_x == pos_end_y
        bk_pt_overlap = "Right"
    else
        bk_pt_overlap = "Ranged"
    end
    consensus = string(split(split(y)[11],":")[10] == split(split(y)[10],":")[10])
    bk_pt_overlap = bk_pt_overlap*"\t"*consensus*"\t"*replace(y,"\t" => "|") 
    return bk_pt_overlap
end

function ann_sv(sv,l1, l2, comp_col)
    if sv =="Unique"
        return l1*"\t"*sv*"\tNA\tNA\tNA"
    else
        overlap = compare_break_points(l1,l2[1],comp_col)
        #println("Done3")
        return l1*"\t"*sv*"\t"*overlap
    end
end

function annotate_pindel(v1,v2,dat_col,comp_col)
    check_line = iterate(v2)
    v3=[]
    #i = 1
    #j = 1
    for line1 in v1
        if check_pos(line1,check_line) 
            sv = comp_loc(line1,check_line[1],dat_col)
        elseif iterate(v2, check_line[2]) !== nothing
            check_line = iterate(v2,check_line[2])
            #j +=1
            sv= comp_loc(line1,check_line[1],dat_col)
        else
            sv= comp_loc(line1,check_line[1],dat_col)
        end
        val = ann_sv(sv,line1,check_line,comp_col)
        #println(val)
        push!(v3,val)
        #i +=1
        #println(string(i)*" " *string(j))
    end
    return v3
end

function print_res(filePath,vec)
    fheader = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT\tAD_REF\tAD_ALT"
    fheader = fheader*"\tStatus_of_merge\tbreak_point_kind\tSeq_overlap_status\tResults_of_merge"
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
        tmpline_to_print[14] = first(tmpline_to_print[14],10)

        push!(tmpline_to_print, find_pos_val("SVTYPE=",line))
        push!(tmpline_to_print, find_pos_val("SVLEN=",line))
        push!(tmpline_to_print, find_pos_val("HOMLEN=",line))
        push!(tmpline_to_print, string(AD_alt/(AD_alt+AD_ref)))
        println(f,join(tmpline_to_print,"\t"))
    end
    close(f)
end

datdir = "/Users/aragaven/scratch/test_julia/test_pindel/pindel_vcf";
samp1 = "DSY303_OP";
samp2 = "DSY303_WH";
merged_file = "/Users/aragaven/Documents/Research/CBC/ene/DSY303_merged.vcf";
vec1 = readlines(joinpath(datdir,samp1*".vcf"));
vec2 = readlines(joinpath(datdir,samp2*".vcf"));
vec3 = readlines(merged_file);
vec1_1 = filter(x -> !startswith(x, "#"), vec1);
vec2_2 = filter(x -> !startswith(x, "#"), vec2);
vec3_3 = filter(x -> !startswith(x, "#"), vec3);
vec4 = annotate_pindel(vec1_1,vec3_3,10,11);
vec5 = annotate_pindel(vec2_2,vec3_3,11,10);

filep1 = "/Users/aragaven/Documents/Research/CBC/ene/DSY303_OP.results.tsv"
print_res(filep1,vec4)
filep2 = "/Users/aragaven/Documents/Research/CBC/ene/DSY303_WH.results.tsv"
print_res(filep2,vec5)

