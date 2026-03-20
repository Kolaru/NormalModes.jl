ORCA_UNITS = Dict(
    "hessian" => aunit(u"J/m^2"),
    "vibrational_frequencies" => u"cm^-1",
    "positions" => u"bohr"
)

function load_orca(file)
    data = Dict()

    for p in split(read(file, String), '$')
        lines = filter(
            line -> (length(line) > 0 && line[1] != '#'),
            split(p, '\n')
        )
        isempty(lines) && continue

        key = lines[1]

        if key in ["act_atom", "act_coord", "act_energy", "multiplicity",
                "actual_temperature", "frequency_scale_factor"]
            data[key] = parse(Float64, strip(lines[2]))

        elseif key in ["hessian", "normal_modes", "vibrational_frequencies", "dipole_derivatives", "ir_spectrum", "atoms"]
            n = parse(Int, first(split(lines[2])))
            matrix = Matrix{String}(undef, n, 0)
            colstart = 1
            multiline = length(lines) > 2n

            for (k, line) in enumerate(lines[3:end])
                row = mod1(k, n + 1)

                if multiline
                    row -= 1
                    row == 0 && continue
                    cells = split(line)[2:end]
                else
                    cells = split(line)[1:end]
                end

                if row == 1
                    matrix = hcat(matrix, Matrix{String}(undef, n, length(cells)))
                end

                matrix[row, colstart:(colstart - 1 + length(cells))] = cells

                if row == n
                    colstart += length(cells)
                end
            end

            if key == "atoms"
                data["atoms"] = Symbol.(matrix[:, 1])
                data["positions"] = parse.(Float64, matrix[:, 2:end])
            elseif key == "vibrational_frequencies"
                data["vibrational_frequencies"] = parse.(Float64, matrix[:, 2])
            else
                matrix = parse.(Float64, matrix)

                if size(matrix, 2) == 1
                    data[key] = vec(matrix)
                else
                    data[key] = matrix
                end
            end
        end
    end

    for (key, unit) in ORCA_UNITS
        data[key] = data[key] * unit
    end

    return data
end

function load_orca(::Type{NormalDecomposition}, file)
    data = load_orca(file)
    return NormalDecomposition(data["hessian"], data["atoms"])
end