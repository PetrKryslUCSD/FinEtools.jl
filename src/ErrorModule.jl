module ErrorModule
using ..NodalFieldModule: NodalField
using ..ElementalFieldModule: ElementalField
using ..FEMMBaseModule: AbstractFEMM
using ..IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using ..FieldModule: gathervalues_asmat!, ndofs
using ..MatrixUtilityModule: locjac!
using ..FESetModule: nodesperelem, manifdim


export L2error

function L2error(self::AbstractFEMM,
                 geom::NodalField{GFT},
                 u::NodalField{UFT}, 
                 exact::Function) where {GFT,UFT}
    
    fes = self.integdomain.fes
    dimu = size(u.values, 2)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)

    fes = self.integdomain.fes
    ndn = ndofs(u) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element

    ecoords = fill(zero(GFT), nne, sdim) # array of Element coordinates
    loc = fill(zero(GFT), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(GFT), sdim, mdim) # Jacobian matrix -- buffer

    err = ElementalField(zeros(count(fes), 1))
    for i in 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        elvec = u.values[[fes.conn[i]...],:]
        for j in 1:npts
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            err.values[i] += ((elvec'*Ns[j] .- exact(ecoords))' * (elvec'*Ns[j] .- exact(ecoords)) * w[j] * Jac)[1]
        end
        err.values[i] = sqrt(err.values[i])
    end
    return err
end


end
