import math
import numpy as np
from stl import mesh

# Input
length = 2.17
definition = 5
yDef = 50
radius = []
for x in range(yDef):
    x = x/yDef*length
    radius.append(0.5 + (1.3*(x-1.06)/2.17)**3 + (0.65*x / 2.17)**2 - 0.3*x)
minRadius = 0.2
cut = True
efficient = True
title = 'low_poly.stl'

# Scaling
scalingFactor = 25.4
length *= scalingFactor
minRadius *= scalingFactor
radius = [x*scalingFactor for x in radius]

# Check input bounds
if min(radius) < minRadius:
    raise ValueError('minRadius cannot be greater than any radii values.')
if definition <= 1:
    raise ValueError('definition cannot be <= 1.')
if length < 0:
    raise ValueError('length cannot be negative.')
if not radius:
    raise ValueError('radius cannot be empty.')
if minRadius < 0:
    raise ValueError('minRadius cannot be negative.')  

# Creates vertices
cutDef = int(definition if not cut else definition/2 + 1)

    # External vertices
externalV = []
for i, r in enumerate(radius):
    for j in range(cutDef):
        x = i / (len(radius) - 1) * length
        y = math.sin(j / definition * 2*math.pi) * r
        z = math.cos(j / definition * 2*math.pi) * r
        externalV.append([x, -z, y])

    # Internal vertices
internalV = []
if minRadius == 0:
    internalV.append([0, 0, 0])
    internalV.append([length, 0, 0])
else:
    if cut:
        for i, r in enumerate(radius):
            for j in range(cutDef):
                x = i / (len(radius) - 1) * length
                y = math.sin(j / definition * 2*math.pi) * minRadius
                z = math.cos(j / definition * 2*math.pi) * minRadius
                internalV.append([x, -z, y])
    else:
        for i in range(cutDef):
            y = math.sin(i / definition * 2*math.pi) * minRadius
            z = math.cos(i / definition * 2*math.pi) * minRadius
            internalV.append([0, -z, y])
        for i in range(cutDef):
            y = math.sin(i / definition * 2*math.pi) * minRadius
            z = math.cos(i / definition * 2*math.pi) * minRadius
            internalV.append([length, -z, y])
    
vertices = np.concatenate((externalV, internalV))



def meshify(bl, br, tr, tl):
    return [[bl, tr, tl], [bl, br, tr]]

# Creates faces
externalF = []
for i in range(len(radius) - 1):
    for j in range(cutDef):
        if not cut or j != cutDef - 1:
            a = j + 1 if j != cutDef - 1 else 0 # Wraps the "seam" of the rotated solid
            externalF.extend(meshify(
                i*cutDef + a,
                i*cutDef + j,
                (i + 1)*cutDef + j,
                (i + 1)*cutDef + a
                ))
       
# Internal and side faces
sideF = []
if minRadius == 0:
    # Side faces
    for i in range(cutDef):
        if not cut or i != cutDef - 1:
            a = i + 1 if i != cutDef - 1 else 0
            sideF.append([
                -2,
                i,
                a
                ])
            sideF.append([
                -1,
                len(externalV) - cutDef + a,
                len(externalV) - cutDef + i
                ])
else:
    # Internal faces
    indexOffset = -cutDef if cut else len(externalV) + cutDef
    internalF = []
    if efficient:
        for i in range(cutDef): # Simpler, less faces, may not work
            if not cut or i != cutDef - 1:
                a = i + 1 if i != cutDef - 1 else 0
                internalF.extend(meshify(
                    len(externalV) + i,
                    len(externalV) + a,
                    indexOffset + a,
                    indexOffset + i
                    ))
    else:
        for i in range(len(radius) - 1): # More faces, likely to work
            for j in range(cutDef):
                if not cut or j != cutDef - 1:
                    a = j + 1 if j != cutDef - 1 else 0
                    internalF.extend(meshify(
                        len(externalV) + i*cutDef + j,
                        len(externalV) + i*cutDef + a,
                        len(externalV) + (i + 1)*cutDef + a,
                        len(externalV) + (i + 1)*cutDef + j
                        ))
    # Side faces
    for i in range(cutDef):
        if not cut or i != cutDef - 1:
            a = i + 1 if i != cutDef - 1 else 0
            sideF.extend(meshify(
                len(externalV) + a,
                len(externalV) + i,
                i,
                a
                ))
            sideF.extend(meshify(
                indexOffset + i,
                indexOffset + a,
                len(externalV) - cutDef + a,
                len(externalV) - cutDef + i
                ))

faces = np.concatenate((externalF, sideF))
if minRadius != 0:
    faces = np.concatenate((internalF, faces))

# Cut faces
if cut:
    cutF = []
    for i in range(len(radius) - 1):
        cutF.extend(meshify(
            len(externalV) + i*cutDef,
            len(externalV) + (i + 1)*cutDef,
            (i + 1)*cutDef,
            i*cutDef
            ))
        cutF.extend(meshify(
            (i + 1)*cutDef - 1,
            (i + 2)*cutDef - 1,
            len(externalV) + (i + 2)*cutDef - 1,
            len(externalV) + (i + 1)*cutDef - 1,
            ))
    faces = np.concatenate((faces, cutF))

# Create the mesh
cube = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        cube.vectors[i][j] = vertices[f[j],:]
cube.save(title)

print(f'Vertices: {len(vertices)}')
print(f'Faces: {len(faces)}')