/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an
    unofficial extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "isoCutFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::isoCutFace, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoCutFace::isoCutFace
(
    const fvMesh& mesh,
    scalarField& f
)
:
    mesh_(mesh),
    f_(f),
    firstEdgeCut_(-1),
    lastEdgeCut_(-1),
    firstFullySubmergedPoint_(-1),
    nFullySubmergedPoints_(0),
    subFaceCentre_(vector::zero),
    subFaceArea_(vector::zero),
    subFacePoints_(10),
    surfacePoints_(4),
    subFaceCentreAndAreaIsCalculated_(false)
{
    clearStorage();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::isoCutFace::calcSubFaceCentreAndArea()
{
    const label nPoints = subFacePoints_.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        subFaceCentre_ = sum(subFacePoints_)/scalar(3);
        subFaceArea_ =
            0.5
           *(
                (subFacePoints_[1] - subFacePoints_[0])
               ^(subFacePoints_[2] - subFacePoints_[0])
            );
    }
    else if (nPoints > 0)
    {
        vector sumN = vector::zero;
        scalar sumA = 0.0;
        vector sumAc = vector::zero;
        const point fCentre = sum(subFacePoints_)/scalar(nPoints);

        for (label pI = 0; pI < nPoints; ++pI)
        {
            const point& nextPoint = subFacePoints_[subFacePoints_.fcIndex(pI)];

            vector c = subFacePoints_[pI] + nextPoint + fCentre;
            vector n =
                (nextPoint - subFacePoints_[pI])^(fCentre - subFacePoints_[pI]);
            scalar a = magSqr(n);

            sumN += n;
            sumA += a;
            sumAc += a*c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA < ROOTVSMALL)
        {
            subFaceCentre_ = fCentre;
            subFaceArea_ = vector::zero;
        }
        else
        {
            subFaceCentre_ = (1.0/3.0)*sumAc/sumA;
            subFaceArea_ = 0.5*sumN;
        }
    }

    subFaceCentreAndAreaIsCalculated_ = true;
}


Foam::label Foam::isoCutFace::calcSubFace
(
    const scalar isoValue,
    const pointField& points,
    const scalarField& f,
    const labelList& pLabels
)
{
    // Face status set to one of the values:
    //  -1: face is fully below the isosurface
    //   0: face is cut, i.e. has values larger and smaller than isoValue
    //  +1: face is fully above the isosurface
    label faceStatus;

    label pl1 = pLabels[0];
    scalar f1 = f[pl1];

    // If vertex values are very close to isoValue lift them slightly to avoid
    // dealing with the many special cases of a face being touched either at a
    // single point, along an edge, or the entire face being on the surface.
    if (mag(f1 - isoValue) < 10*SMALL)
    {
        f1 += sign(f1 - isoValue)*10*SMALL;
    }

    // Finding cut edges, the point along them where they are cut, and all fully
    // submerged face points.
    forAll(pLabels, pI)
    {
        label pl2 = pLabels[pLabels.fcIndex(pI)];
        scalar f2 = f[pl2];
        if (mag(f2 - isoValue) < 10*SMALL)
        {
            f2 += sign(f2 - isoValue)*10*SMALL;
        }

        if (f1 > isoValue)
        {
            nFullySubmergedPoints_ += 1;

            if (f2 < isoValue)
            {
                lastEdgeCut_ = (isoValue - f1)/(f2 - f1);
            }
        }
        else if (f1 < isoValue && f2 > isoValue)
        {
            if (firstFullySubmergedPoint_ == -1)
            {
                firstFullySubmergedPoint_ = pLabels.fcIndex(pI);

                firstEdgeCut_ = (isoValue - f1)/(f2 - f1);
            }
            else
            {
                if (debug)
                {
                    const face fl(pLabels);

                    WarningIn
                    (
                         "Foam::label Foam::isoCutFace::calcSubFace(...)"
                    )   << "More than two face cuts for face " << fl
                        << endl;

                    Pout<< "Face values: f-isoValue = " << endl;
                    forAll(fl, fpi)
                    {
                        Pout<< f[fl[fpi]] - isoValue << " ";
                    }
                    Pout<< " " << endl;
                }
            }
        }
        pl1 = pl2;
        f1 = f2;
    }

    if (firstFullySubmergedPoint_ != -1)
    {
        // Face is cut
        faceStatus = 0;
        subFacePoints(points, pLabels);
    }
    else if (f1 < isoValue)
    {
        // Face entirely above isosurface
        faceStatus = 1;
    }
    else
    {
        // Face entirely below isosurface
        faceStatus = -1;
    }

    return faceStatus;
}


void Foam::isoCutFace::subFacePoints
(
    const pointField& points,
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();

    surfacePoints(points, pLabels);

    forAll(surfacePoints_, pI)
    {
        subFacePoints_.append(surfacePoints_[pI]);
    }

    for (label pI = 0; pI < nFullySubmergedPoints_; ++pI)
    {
        subFacePoints_.append
        (
            points[pLabels[(firstFullySubmergedPoint_ + pI) % nPoints]]
        );
    }
}


void Foam::isoCutFace::surfacePoints
(
    const pointField& points,
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();

    const label n = firstFullySubmergedPoint_ + nFullySubmergedPoints_;

    label pl1 = pLabels[(n - 1) % nPoints];

    label pl2 = pLabels[n % nPoints];

    surfacePoints_.append
    (
        points[pl1] + lastEdgeCut_*(points[pl2] - points[pl1])
    );

    pl1 = pLabels[(firstFullySubmergedPoint_ - 1 + nPoints) % nPoints];
    pl2 = pLabels[firstFullySubmergedPoint_];

    surfacePoints_.append
    (
        points[pl1] + firstEdgeCut_*(points[pl2] - points[pl1])
    );
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::isoCutFace::calcSubFace
(
    const label faceI,
    const scalar isoValue
)
{
    clearStorage();
    const labelList& pLabels = mesh_.faces()[faceI];
    const pointField& points = mesh_.points();
    return calcSubFace(isoValue, points, f_, pLabels);
}


Foam::label Foam::isoCutFace::calcSubFace
(
    const pointField& points,
    const scalarField& f,
    const scalar isoValue
)
{
    clearStorage();
    const labelList pLabels(identity(f.size()));
    return calcSubFace(isoValue, points, f, pLabels);
}


const Foam::point& Foam::isoCutFace::subFaceCentre()
{
    if (!subFaceCentreAndAreaIsCalculated_)
    {
        calcSubFaceCentreAndArea();
    }
    return subFaceCentre_;
}


const Foam::vector& Foam::isoCutFace::subFaceArea()
{
    if (!subFaceCentreAndAreaIsCalculated_)
    {
        calcSubFaceCentreAndArea();
    }
    return subFaceArea_;
}


const Foam::DynamicList<Foam::point>& Foam::isoCutFace::subFacePoints() const
{
    return subFacePoints_;
}


const Foam::DynamicList<Foam::point>& Foam::isoCutFace::surfacePoints() const
{
    return surfacePoints_;
}


void Foam::isoCutFace::clearStorage()
{
    firstEdgeCut_ = -1;
    lastEdgeCut_ = -1;
    firstFullySubmergedPoint_ = -1;
    nFullySubmergedPoints_ = 0;
    subFaceCentre_ = vector::zero;
    subFaceArea_ = vector::zero;
    subFacePoints_.clear();
    surfacePoints_.clear();
    subFaceCentreAndAreaIsCalculated_ = false;
}


Foam::scalar Foam::isoCutFace::timeIntegratedArea
(
    const pointField& fPts,
    const scalarField& pTimes,
    const scalar dt,
    const scalar magSf,
    const scalar Un0
)
{
    // Initialise time integrated area returned by this function
    scalar tIntArea = 0.0;

    // Finding ordering of vertex points
    labelList order(pTimes.size());
    sortedOrder(pTimes, order);
    const scalar firstTime = pTimes[order.first()];
    const scalar lastTime = pTimes[order.last()];

    // Dealing with case where face is not cut by surface during time interval
    // [0,dt] because face was already passed by surface
    if (lastTime <= 0)
    {
        // If all face cuttings were in the past and cell is filling up (Un0>0)
        // then face must be full during whole time interval
        tIntArea = magSf*dt*pos(Un0);
        return tIntArea;
    }

    // Dealing with case where face is not cut by surface during time interval
    // [0, dt] because dt is too small for surface to reach closest face point
    if (firstTime >= dt)
    {
        // If all cuttings are in the future but non of them within [0,dt] then
        // if cell is filling up (Un0 > 0) face must be empty during whole time
        // interval
        tIntArea = magSf*dt*(1 - pos(Un0));
        return tIntArea;
    }

    // If we reach this point in the code some part of the face will be swept
    // during [tSmall, dt-tSmall]. However, it may be the case that there are no
    // vertex times within the interval. This will happen sometimes for small
    // time steps where both the initial and the final face-interface
    // intersection line (FIIL) will be along the same two edges.

    // Face-interface intersection line (FIIL) to be swept across face
    DynamicList<point> FIIL(3);
    // Submerged area at beginning of each sub time interval time
    scalar initialArea = 0.0;
    //Running time keeper variable for the integration process
    scalar time = 0.0;

    // Special treatment of first sub time interval
    if (firstTime > 0)
    {
        // If firstTime > 0 the face is uncut in the time interval
        // [0, firstTime] and hence fully submerged in fluid A or B.
        // If Un0 > 0 cell is filling up and it must initially be empty.
        // If Un0 < 0 cell must initially be full(y immersed in fluid A).
        time = firstTime;
        initialArea = magSf*(1.0 - pos(Un0));
        tIntArea = initialArea*time;
        cutPoints(fPts, pTimes, time, FIIL);
    }
    else
    {
        // If firstTime <= 0 then face is initially cut and we must
        // calculate the initial submerged area and FIIL:
        time = 0.0;
        // Note: calcSubFace assumes well-defined 2-point FIIL!!!!
        calcSubFace(fPts, -sign(Un0)*pTimes, time);
        initialArea = mag(subFaceArea());
        cutPoints(fPts, pTimes, time, FIIL);
    }

    // Making sorted array of all vertex times that are between max(0,firstTime)
    // and dt and further than tSmall from the previous time.
    DynamicList<scalar> sortedTimes(pTimes.size());
    {
        scalar prevTime = time;
        const scalar tSmall = max(1e-6*dt, 10*SMALL);
        forAll(order, ti)
        {
            const scalar timeI = pTimes[order[ti]];
            if ( timeI > prevTime + tSmall && timeI <= dt)
            {
                sortedTimes.append(timeI);
                prevTime = timeI;
            }
        }
    }

    // Sweeping all quadrilaterals corresponding to the intervals defined above
    forAll(sortedTimes, ti)
    {
        const scalar newTime = sortedTimes[ti];
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints(fPts, pTimes, newTime, newFIIL);

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea += (newTime - time)*
            (initialArea + sign(Un0)*(alpha/3.0 + 0.5*beta));
        // Adding quad area to submerged area
        initialArea += sign(Un0)*(alpha + beta);

        FIIL = newFIIL;
        time = newTime;
    }

    // Special treatment of last time interval
    if (lastTime > dt)
    {
        // FIIL will end up cutting the face at dt
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints(fPts, pTimes, dt, newFIIL);

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea += (dt - time)*
            (initialArea + sign(Un0)*(alpha/3.0 + 0.5*beta));
    }
    else
    {
        // FIIL will leave the face at lastTime and face will be fully in fluid
        // A or fluid B in the time interval from lastTime to dt.
        tIntArea += magSf*(dt - lastTime)*pos(Un0);
    }

    return tIntArea;
}


void Foam::isoCutFace::cutPoints
(
    const pointField& pts,
    const scalarField& f,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    const label nPoints = pts.size();
    scalar f1(f[0]);

    // Snapping vertex value to f0 if very close (needed for 2D cases)
    if (mag(f1 - f0) < 10*SMALL)
    {
        f1 = f0;
    }

    forAll(pts, pI)
    {
        label pI2 = (pI + 1) % nPoints;
        scalar f2 = f[pI2];

        // Snapping vertex value
        if (mag(f2 - f0) < 10*SMALL)
        {
            f2 = f0;
        }

        if ((f1 < f0 && f2 > f0) || (f1 > f0 && f2 < f0))
        {
            const scalar s = (f0 - f1)/(f2 - f1);
            cutPoints.append(pts[pI] + s*(pts[pI2] - pts[pI]));
        }
        else if (f1 == f0)
        {
            cutPoints.append(pts[pI]);
        }
        f1 = f2;
    }

    if (cutPoints.size() > 2)
    {
        WarningIn
        (
            "void Foam::isoCutFace::cutPoints(...)"
        )   << "cutPoints = " << cutPoints << " for pts = " << pts
            << ", f - f0 = " << f - f0 << " and f0 = " << f0 << endl;
    }
}


void Foam::isoCutFace::quadAreaCoeffs
(
    const DynamicList<point>& pf0,
    const DynamicList<point>& pf1,
    scalar& alpha,
    scalar& beta
) const
{
    // Number of points in provided face-interface intersection lines
    const label np0 = pf0.size();
    const label np1 = pf1.size();

    // quad area coeffs such that area(t) = alpha*t^2 + beta*t.
    // With time interval normalised, we have full quadArea = alpha + beta
    // and time integrated quad area = alpha/3 + beta/2;
    alpha = 0.0;
    beta = 0.0;

    if (np0 && np1)
    {
        // Initialising quadrilateral vertices A, B, C and D
        vector A(pf0[0]);
        vector C(pf1[0]);
        vector B(pf0[0]);
        vector D(pf1[0]);

        if (np0 == 2)
        {
            B = pf0[1];
        }
        else if (np0 > 2)
        {
            WarningIn
            (
                "void Foam::isoCutFace::quadAreaCoeffs(...)"
            )   << "Vertex face was cut at pf0 = " << pf0 << endl;
        }

        if (np1 == 2)
        {
            D = pf1[1];
        }
        else if (np1 > 2)
        {
            WarningIn
            (
                "void Foam::isoCutFace::quadAreaCoeffs(...)"
            )    << "Vertex face was cut at pf1 = " << pf1 << endl;
        }

        // Swapping pf1 points if pf0 and pf1 point in same general direction
        // (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if (((B - A) & (D - C)) > 0)
        {
            vector tmp = D;
            D = C;
            C = tmp;
        }

        // Defining local coordinates (xhat, yhat) for area integration of swept
        // quadrilateral ABCD such that A = (0,0), B = (Bx,0), C = (Cx,Cy) and
        // D = (Dx,Dy) with Cy = 0 and Dy > 0.

        const scalar Bx = mag(B - A);

        vector xhat(vector::zero);
        if (Bx > 10*SMALL)
        {
            // If |AB| > 0 ABCD we use AB to define xhat
            xhat = (B - A)/mag(B - A);
        }
        else if (mag(C - D) > 10*SMALL)
        {
            // If |AB| ~ 0 ABCD is a triangle ACD and we use CD for xhat
            xhat = (C - D)/mag(C - D);
        }
        else
        {
            return;
        }

        // Defining vertical axis in local coordinates
        vector yhat = D - A;
        yhat -= ((yhat & xhat)*xhat);

        if (mag(yhat) > 10*SMALL)
        {
            yhat /= mag(yhat);

            const scalar Cx = (C - A) & xhat;
            const scalar Cy = mag((C - A) & yhat);
            const scalar Dx = (D - A) & xhat;
            const scalar Dy = mag((D - A) & yhat);

            // area = ((Cx - Bx)*Dy - Dx*Cy)/6.0 + 0.25*Bx*(Dy + Cy);
            alpha = 0.5*((Cx - Bx)*Dy - Dx*Cy);
            beta = 0.5*Bx*(Dy + Cy);
        }
    }
    else
    {
        WarningIn
        (
            "void Foam::isoCutFace::quadAreaCoeffs(...)"
        )   << "Vertex face was cut at " << pf0 << " and at " << pf1 << endl;
    }
}


// ************************************************************************* //
