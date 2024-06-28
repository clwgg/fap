using Pkg
# Pkg.add("DataFrames")
# Pkg.add("Distributions")

using DataFrames
using Random
using Distributions

function mutationFrequency(cellGenomeFull, pop_size, seed_origins, raw_threshold=0.01)
    # Flatten mutations
    muts = vcat(cellGenomeFull...)

    # Count mutations
    mutation_counts = countmap(muts)

    # Mutation frequency
    mutation_frequency = Dict{Int, Float64}()
    for (mutation, count) in pairs(mutation_counts)
        mutation_frequency[mutation] = count / (pop_size * 2.0) # diploid
    end

    # Remove mutations with frequency less than threshold
    mutation_frequency = Dict(filter(x -> x[2] >= raw_threshold, mutation_frequency))

    # Rebuild seed_origins dict
    seed_origins_return = Dict{Int, Int}()
    # loop over mutations in the frequency dict
    for (mutation, frequency) in mutation_frequency
        seed_origins_return[mutation] = get(seed_origins, mutation, -1)
    end

    # println("Counts complete.")
    return mutation_frequency, seed_origins_return
end

# function mutationFrequency(cellGenomeFull, pop_size, raw_threshold=0.01)
#     # Flatten mutations
#     muts = vcat(cellGenomeFull...)

#     # Unique mutations
#     unique_muts = unique(muts)

#     # count mutations
#     mutation_counts = Dict{Int, Int}()
#     for mut in unique_muts
#         count = count_instances(mut, cellGenomeFull)
#         # println("Mutation: ", mut, " Count: ", count)
#         mutation_counts[mut] = count
#     end

#     # mutation_counts = Dict{Int, Int}()
#     # for genome in cellGenomeFull
#     #     for mutation in genome
#     #         mutation_counts[mutation] = get(mutation_counts, mutation, 0) + 1
#     #     end
#     # end

#     mutation_frequency = Dict{Int, Float64}()
#     for (mutation, count) in mutation_counts
#         mutation_frequency[mutation] = count / (pop_size * 2.) # diploid
#     end

#     # Remove all mutations with frequency less than threshold
#     mutation_frequency = Dict(filter(x -> x[2] >= raw_threshold, mutation_frequency))

#     println("Counts complete.")
#     return mutation_frequency
# end

function count_instances(value, list_of_lists)
    count = 0
    for inner_list in list_of_lists
        count += sum(x -> x == value, inner_list)
    end
    return count
end

function build_output_table(cellGenomes, mutInductionTimes, driver_muts, group, seed_origin, mean_depth, sd_depth, sample_size=15000)
    # Get population size
    pop_size = length(cellGenomes)

    # Sample the population
    if sample_size != 0
        sample_size = min(sample_size, pop_size)
        sample_idx = sample(1:pop_size, sample_size, replace=false)
        cellGenomes = cellGenomes[sample_idx]
        pop_size = length(cellGenomes)
    end
    
    raw_vafs, seed_origins = mutationFrequency(cellGenomes, pop_size, seed_origin)
    
    output_table = DataFrame(mutation = Int[], raw_vaf = Float64[], vaf = Float64[], t= Int[], is_driver = Int[], seed_parent = Int[], group = String[])
    
    for (mutation, frequency) in raw_vafs
        t = get(mutInductionTimes, mutation, -1) # -1 if can't find the mutation induction time in dict
        is_driver = 0
        for d in 1:length(driver_muts)
            if driver_muts[d] == mutation
                is_driver = 1
            else
                is_driver = 0
                # println("Driver mutation found: ", mutation)
            end
        end

        # Get seq_vafs
        vaf = sequence_mutations(frequency, mean_depth, sd_depth)
        seed_num = get(seed_origins, mutation, -1)

        push!(output_table, (mutation, frequency, vaf, t, is_driver, seed_num, group))
    end
    return output_table
end

function sequence_mutations(vaf, mean_depth, sd_depth)
    # Simulate sequencing mutations getting reads, and depth
    
    # Reads from a binomial distribution with vaf as probability and depth as number of trials
    # Depth is normal distribution with mean of 'depth' and standard deviation 10
    # Return the number of reads
    total_reads = Int(floor(rand(Normal(mean_depth, sd_depth))))

    if total_reads < 0 || vaf > 1 || vaf < 0
        println("Total Reads: ", total_reads)
        println("VAF: ", vaf)
        println("ERROR!")
        exit()
    end

    var_reads_dist = Binomial( total_reads, vaf )
    var_reads = rand(var_reads_dist)

    # var_reads = rand(var_reads_dist, rand_total_reads)

    # Divide each var_reads by the total reads
    return var_reads / total_reads
end

function assign_seed_origin(polyp_cellGenomes::Vector{Vector{Int}}, seed_genomes::Vector{Vector{Int}})
    # Create a dictionary to store the origin of each mutation
    origindict = Dict{Int, Int}()

    # Convert seed genomes to sets for faster operations
    seed_sets = [Set(seed) for seed in seed_genomes]

    # Loop through each cell genome
    for (i, cell_genome) in enumerate(polyp_cellGenomes)
        cell_set = Set(cell_genome)
        # Loop through each seed genome set
        for (j, seed_set) in enumerate(seed_sets)
            # Check if the seed genome is a subset of the cell genome
            if issubset(seed_set, cell_set)
                # Assign the origin for each mutation in the cell genome
                for m in cell_genome
                    origindict[m] = j
                end
                break # Once matched, no need to check other seed genomes
            end
        end
    end

    return origindict
end