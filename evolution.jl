type Individual
    valid_genome::Array{Float64, 1}  # genome used for validation
    chromosome::Array{Float64, 1}
    mutation_chance::Float64
    function Individual(valid_genome, chromosome::Array{Float64, 1}, mutation_chance=0.1)
        this = new(valid_genome, chromosome, mutation_chance)
        repair!(this)
        return this
    end
    

end

function mutate!(self::Individual)
    if rand() >= self.m_chance
        return
    end
    i = rand(1:length(self.chromosome) - 3)
    j = rand(i+1:length(self.chromosome) - 2)
    center = self.chromosome[i:j]
    reverse!(center)
    self.chromosome = vcat(self.chromosome[1:i], center, self.chromosome[j:end])
end

function mutate(self::Individual)
    chromosome_copy = self.chromosome
    mutate!(self)
    another = Individual(self.valid_genome, self.chromosome, self.mutation_chance)
    self.chromosome = chromosome_copy;
    return another
end

function crossover(self::Individual, other::Individual)
    chromosome1, chromosome2 = _cross_over_chromosomes(self, other)
    i1 = Individual(self.valid_genome, chromosome1, self.m_chance)
    i2 = Individual(self.valid_genome, chromosome2, self.m_chance)
    return i1, i2
end

function _cross_over_chromosomes(ch1::Array{Float64}, ch2::Array{Float64}) #krzyżowanie chromosomów 29216.4532612964 27748.70957813485
    i = rand(1:length(ch1) - 3)
    j = rand(i+1:length(ch1) - 2)
    # dla pierwszego dziecka
    middle = ch1[i:j]  # pobierz kawalek genomu rodzica 1
    tmp = filter(x -> !(x in middle), ch2)  # wyrzuc z genomu 2 geny wystepujace w middle

    start, ending = tmp[1:i-1], tmp[i:end]  # podziel pozostalosc genomu 2 na 2 czesci pasujace do genomu 1
    child1 = vcat(start, middle, ending) # polacz czesci razem
    #dla drugiego dziecka
    middle = ch2[i:j]
    tmp = filter(x -> !(x in middle), ch1)
    start, ending = tmp[1:i-1], tmp[i:end]
    child2 = vcat(start, middle, ending)
    return child1, child2
end

function repair!(self::Individual)  # gets rid of duplicates
    doubled = Dict(Float64, Int)
    for elem in self.chromosome  # build dict with <elem, number of occuriences> pairs 
        if haskey(doubled, elem)
            doubled[elem] += 1
        else
            doubled[elem] = 1
        end
    end
    filter!((key, val) -> val > 1, doubled)
    
    keys = vcat(keys(doubled))
    chromosome_set = Set(self.chromosome)
    for elem in self.valid_genome
        if elem in chromosome_set  # that elem is present with genome of this individual, search for next one
            continue
        end
        self.chromosome[getindex(self.chromosome, keys[end])] = i  #podmień i za powtarzający się genom        
    
        doubled[keys[end]] -= 1
        if doubled[keys[end]] == 0
            pop!(doubled, keys[end])
            deleteat!(keys, keys[end])
        end
    end
end
print("wat")