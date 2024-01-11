using DataFrames, CSV, StatsBase, Base.Threads

function read_plink(file)
    bed_file = string(file, ".bed")
    bim_file = string(file, ".bim")
    fam_file = string(file, ".fam")
    bim = CSV.read(bim_file, DataFrame; header = false)
    rename!(bim, [:Chromosome, :SNP_ID, :Genetic_Distance, :Physical_Position, :Allele1, :Allele2])
    fam = CSV.read(fam_file, DataFrame; header = false, missingstring="NA")
    # rename!(fam[:, 1:5], [:Family_ID, :Individual_ID, :Paternal_ID, :Maternal_ID, :Sex])
    num_phenotypes = size(fam, 2) - 5  # 表型的数量
    traits = Symbol.("Trait_" .* string.(1:num_phenotypes))
    rename!(fam, [:Family_ID, :Individual_ID, :Paternal_ID, :Maternal_ID, :Sex, traits...])
    bed = read_bed(bed_file, size(fam, 1), size(bim, 1))
    return (bim = bim, fam = fam, bed = bed)
    # return bim
end


function read_bed(genfil, nind, nloci)

    code = Dict(0 => 2, 2 => 1, 1 => 3, 3 => 0) # an array for mapping genotype values

    # open a binary file, read plink magic numbers and mode, check if it is SNP-major mode
    f = open(genfil, "r")
    b1 = read(f, UInt8) # plink magic number 1
    if b1 != parse(UInt8,"01101100",base =2)
        println("Binary genotype file may not be a plink file")
        return
    end

    b1 = read(f, UInt8) # plink magic number 2
    if b1 != parse(UInt8,"00011011",base =2)
        println("Binary genotype file may not be a plink file")
        return
    end

    b1 = read(f, UInt8) # mode
    if b1 != 0x01
        println("SNP file should not be in individual mode")
        return
    end

    len = filesize(f) - 3
    buffer = Vector{Int8}(undef, len*4)

    for i in 1:len
        p = read(f, UInt8)
        # r = i ÷ nind
        for x in 0:3
           buffer[4 * (i - 1) + x + 1] = code[(p >> (2 * x))&0x03]
        end
    end

    close(f)
    geno = reshape(buffer, Int(len*4/nloci), nloci)[1:nind, :]
    # return the geno matrix
    return geno
end

# read_bed("D:/Datasets/kaml_example/Plink_binary/testdata.bed", nind, nloci)
# data = MIBLUP.read_plink("D:/Datasets/kaml_example/Plink_binary/testdata")
# data.bed
# mean(data.bed, dims = 1) ./ 2

function geno_imputation(x::Vector; missing_value = 3)
    if in(missing_value, x)
        index = x .!= missing_value
        missing_index = .!index
        vec = x[index]
        max_val = mode(vec)
        x[missing_index] .= max_val
    end
    return x
end

function geno_imputation(X::Matrix; missing_value = 3)
    res = copy(X)
    for i in 1:size(X, 2)
        res[:, i] = geno_imputation(X[:, i]; missing_value = missing_value)
    end
    return res
end

# data.bed
# MIBLUP.geno_imputation(data.bed)
# data.bed
function read_geno_file(filename)
    # 打开文件并读取所有行
    lines = readlines(filename)
    # 去掉每行的首尾空白并按空白分割
    lines = [split(strip(line)) for line in lines]
    # 从第二列开始提取基因型数据，同时保存第一列的值作为行名
    genotypes = []
    rownames = []
    for line in lines
        # 解析基因型数据并验证每个元素是否为整数
        genotype = []
        for x in line[2]
            try
                push!(genotype, parse(Int, x))
            catch e
                error("Invalid genotype data: $x")
            end
        end
        push!(genotypes, genotype)
        push!(rownames, line[1])
    end
    # 验证所有样本的基因型数据长度是否一致
    geno_lengths = Set(length.(genotypes))
    if length(geno_lengths) > 1
        error("Inconsistent genotype lengths")
    end
    # 将基因型数据转换为矩阵并进行转置
    return Matrix(hcat(genotypes...)'), rownames
end
