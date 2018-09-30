//Standard C++:
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
//My Headers:
#include "vector2.h"
#include "materialDatabase.h"

class Particle
{
    private:
        //Coords:
        vector2 coords;
        //Velocities:
        vector2 velocity;
        //Material:
        std::string material;
        //Expanding:
        bool fillToolExpands = false;
        //Mass:
        double mass = 0;
        //Bounciness:
        double bounciness = 0.5;
    public:
        //All values:
        void setAllValues(vector2, vector2, std::string, double, double);
        void setEmpty();

        //Copy Particle:
        void copyParticle(Particle&);

        //Coords:
        void setCoords(vector2);
        vector2 getCoords();
        vector2 getPreciseCoords();

        //Velocities:
        void giveVelocity(vector2);
        void setVelocity(vector2);
        vector2 getVelocity();

        //Material:
        void setMaterial(std::string);
        void setMaterialProperties(std::string);
        std::string getMaterial();

        //Expanding:
        void setFillParticle(bool);
        bool isFillParticle();

        //Mass:
        void setMass(double);
        double getMass();

        //Bounciness:
        void setBounciness(double);
        double getBounciness();

        //Gravitational Velocity:
        void calculateGravitationalVelocity(Particle&);

        //Update:
        void update();
};

//Set values:
void Particle::setAllValues(vector2 startCoords, vector2 startVelocity, std::string startMaterial, double startMass, double startBounciness)
{
    coords = startCoords;
    velocity = startVelocity;
    material = startMaterial;
    mass = startMass;
    bounciness = startBounciness;
}
void Particle::setEmpty()
{
    coords = vector2(floor(coords.x), floor(coords.y));
    velocity = vector2(0, 0);
    material = "empty";
    mass = 0;
    bounciness = 0;
}

//Copy Particle:
void Particle::copyParticle(Particle& particleToCopyTo)
{
    particleToCopyTo.setAllValues(coords, velocity, material, mass, bounciness);
}

//Coords:
void Particle::setCoords(vector2 newCoordinates)
{
    coords = newCoordinates;
}
vector2 Particle::getCoords()
{
    vector2 flooredCoords(floor(coords.x), floor(coords.y));
    return flooredCoords;
}
vector2 Particle::getPreciseCoords()
{
    return coords;
}

//Velocities:
void Particle::giveVelocity(vector2 addedVelocity)
{
    velocity.x = velocity.x + addedVelocity.x;
    velocity.y = velocity.y + addedVelocity.y;
}
void Particle::setVelocity(vector2 newVelocity)
{
    velocity = newVelocity;
}
vector2 Particle::getVelocity()
{
    return velocity;
}

//Material:
void Particle::setMaterial(std::string newMaterial)
{
    material = newMaterial;
}
void Particle::setMaterialProperties(std::string newMaterial)
{
    material = newMaterial;
    mass = getMaterialMass(newMaterial);
    bounciness = getMaterialBounciness(newMaterial);
}
std::string Particle::getMaterial()
{
    return material;
}

//Expanding:
void Particle::setFillParticle(bool isFill)
{
    fillToolExpands = isFill;
}
bool Particle::isFillParticle()
{
    return fillToolExpands;
}

//Mass:
void Particle::setMass(double newMass)
{
    mass = newMass;
}
double Particle::getMass()
{
    return mass;
}

//Bounciness:
void Particle::setBounciness(double newBounciness)
{
    bounciness = newBounciness;
}
double Particle::getBounciness()
{
    return bounciness;
}

//Gravitational Velocity:
void Particle::calculateGravitationalVelocity(Particle& distantParticle)
{
    //Physics constants:
    const double G = 0.00000000006673; //Gravitational Constant (or Big G)

    //Get coords of particle:
    vector2 coords1 = coords;
    //Get coords of particle with gravity:
    vector2 coords2 = distantParticle.getCoords();
    //Get the difference vector:
    vector2 rV(coords2.x - coords1.x, coords2.y - coords1.y);

    //Distances:
    double r = pow(rV.x, 2) + pow(rV.y, 2);
    double r2 = sqrt(r);

    if (r != 0)
    {
        //Normalize the difference vector
        vector2 u(rV.x / r, rV.y / r);
        //Acceleration of gravity
        double a = G * distantParticle.getMass() / r2;
        //Set the velocity:
        velocity.x = velocity.x + (a * u.x / 1000);
        velocity.y = velocity.y + (a * u.y / 1000);
    }
}

//Update:
void Particle::update()
{
    coords.x = coords.x + velocity.x;
    coords.y = coords.y + velocity.y;
}

//Miscellaneous Functions:
void checkParticleMovement(std::vector< std::vector<Particle> >& particleArray)
{
    int vectorWidth = particleArray[0].size();
    int vectorHeight = particleArray.size();
    std::vector< std::vector<bool> > updated(vectorHeight, std::vector<bool> (vectorWidth, 0));

    //Make incrementer:
    int incrementX = 0;
    int incrementY = 0;

    while (incrementY != vectorHeight)
    {
        //Check if it needs to be moved:
        if ((particleArray[incrementY][incrementX].getMaterial() != "empty") && (updated[incrementY][incrementX] == false))
        {
            int coordX = particleArray[incrementY][incrementX].getCoords().x;
            int coordY = particleArray[incrementY][incrementX].getCoords().y;
            //Moving a particle in the grid:
            if ((coordX != incrementX) || (coordY != incrementY))
            {
                if (particleArray[coordY][coordX].getMaterial() == "empty")
                {
                    //Copy Particle:
                    particleArray[incrementY][incrementX].copyParticle(particleArray[coordY][coordX]);
                    //particleArray[coordY][coordX].setCoords(vector2(coordX, coordY));
                    //Delete previous particle:
                    particleArray[incrementY][incrementX].setEmpty();
                    particleArray[incrementY][incrementX].setCoords(vector2(incrementX, incrementY));
                }
            }
            //Make sure the particle can't be updated multiple times:
            updated[coordY][coordX] = true;
        }

        ++incrementX;
        if (incrementX == vectorWidth)
        {
            incrementX = 0;
            ++incrementY;
        }
    }
}

//Collision Detection:
void handleCollisionDetection(std::vector< std::vector<Particle> >& particleArray)
{
    int vectorWidth = particleArray[0].size();
    int vectorHeight = particleArray.size();

    double highestVelocity = 0;

    std::vector< std::vector<vector2> > velocities(vectorHeight, std::vector<vector2>(vectorWidth));
    std::vector< std::vector<vector2> > coords(vectorHeight, std::vector<vector2>(vectorWidth));
    std::vector< std::vector<std::string> > materials(vectorHeight, std::vector<std::string>(vectorWidth));

    //FIND THE HIGHEST VELOCITY (TO DIVIDE WITH):
    int incrementX = 0;
    int incrementY = 0;
    while (incrementY != vectorHeight)
    {
        velocities[incrementY][incrementX] = particleArray[incrementY][incrementX].getVelocity();
        if (velocities[incrementY][incrementX].x > highestVelocity) {highestVelocity = ceil(velocities[incrementY][incrementX].x);}
        if (velocities[incrementY][incrementX].y > highestVelocity) {highestVelocity = ceil(velocities[incrementY][incrementX].y);}

        coords[incrementY][incrementX] = particleArray[incrementY][incrementX].getPreciseCoords();

        materials[incrementY][incrementX] = particleArray[incrementY][incrementX].getMaterial();

        ++incrementX;
        if (incrementX == vectorWidth)
        {
            incrementX = 0;
            ++incrementY;
        }
    }

    //Remove minus number
    highestVelocity = fabs(highestVelocity);

    incrementX = 0;
    incrementY = 0;
    while (incrementY != vectorHeight)
    {
        if (materials[incrementY][incrementX] != "empty")
        {
            vector2 dividedVelocityStart = velocities[incrementY][incrementX];
            if (velocities[incrementY][incrementX].x != 0)
            {
                dividedVelocityStart.x = dividedVelocityStart.x / highestVelocity;
                if (std::isnan(dividedVelocityStart.x) == true) {dividedVelocityStart.x = 0;}
            }
            else {dividedVelocityStart.x = 0;}

            if (velocities[incrementY][incrementX].y != 0)
            {
                dividedVelocityStart.y = dividedVelocityStart.y / highestVelocity;
                if (std::isnan(dividedVelocityStart.y) == true) {dividedVelocityStart.y = 0;}
            }
            else {dividedVelocityStart.y = 0;}

            vector2 dividedVelocityIncrement = dividedVelocityStart;
            while (dividedVelocityIncrement <= velocities[incrementY][incrementX])
            {
                int incrementXLowLimit = incrementX - (highestVelocity * 2);
                if (incrementXLowLimit < 0) {incrementXLowLimit = 0;}
                int incrementXHighLimit = incrementX + (highestVelocity * 2);
                if (incrementXHighLimit >= vectorWidth) {incrementXHighLimit = vectorWidth;}

                int incrementX2 = incrementXLowLimit;

                int incrementY2 = incrementY - (highestVelocity * 2);
                if (incrementY2 < 0) {incrementY2 = 0;}
                int incrementYHighLimit = incrementY + (highestVelocity * 2);
                if (incrementYHighLimit >= vectorHeight) {incrementYHighLimit = vectorHeight;}

                while (incrementY2 != incrementYHighLimit)
                {
                    if ((materials[incrementY2][incrementX2] != "empty") && (incrementX != incrementX2) && (incrementY != incrementY2))
                    {
                        vector2 dividedVelocityStart2 = velocities[incrementY2][incrementX2];
                        vector2 dividedVelocityIncrement2 = dividedVelocityStart2;
                        dividedVelocityIncrement2.x = dividedVelocityIncrement2.x / highestVelocity;
                        if (std::isnan(dividedVelocityIncrement2.x) == true) {dividedVelocityIncrement2.x = 0;}
                        dividedVelocityIncrement2.y = dividedVelocityIncrement2.y / highestVelocity;
                        if (std::isnan(dividedVelocityIncrement2.y) == true) {dividedVelocityIncrement2.y = 0;}

                        while (dividedVelocityIncrement2 <= velocities[incrementY2][incrementX2])
                        {
                            if ((floor(coords[incrementY][incrementX].x + dividedVelocityIncrement.x) ==
                                 floor(coords[incrementY2][incrementX2].x + dividedVelocityIncrement2.x))
                            &&  (floor(coords[incrementY][incrementX].y + dividedVelocityIncrement.y) ==
                                 floor(coords[incrementY2][incrementX2].y + dividedVelocityIncrement2.y)))
                            {
                                std::cout << "COLLISION!" << std::endl;
                            }

                            if (dividedVelocityIncrement2.x >= 0) {dividedVelocityIncrement2.x = dividedVelocityIncrement2.x + dividedVelocityStart2.x;}
                            else {dividedVelocityIncrement2.x = dividedVelocityIncrement2.x - dividedVelocityStart2.x;}

                            if (dividedVelocityIncrement2.y >= 0) {dividedVelocityIncrement2.y = dividedVelocityIncrement2.y + dividedVelocityStart2.y;}
                            else {dividedVelocityIncrement2.y = dividedVelocityIncrement2.y - dividedVelocityStart2.y;}

                            //For minus values:
                            if (dividedVelocityIncrement2 <= 0)
                            {
                                if (dividedVelocityIncrement2 >= velocities[incrementY2][incrementX2]) {break;}
                            }
                        }
                    }

                    ++incrementX2;
                    if (incrementX2 == incrementXHighLimit)
                    {
                        incrementX2 = incrementXLowLimit;
                        ++incrementY2;
                    }
                }

                if (dividedVelocityIncrement.x >= 0) {dividedVelocityIncrement.x = dividedVelocityIncrement.x + dividedVelocityStart.x;}
                else {dividedVelocityIncrement.x = dividedVelocityIncrement.x - dividedVelocityStart.x;}

                if (dividedVelocityIncrement.y >= 0) {dividedVelocityIncrement.y = dividedVelocityIncrement.y + dividedVelocityStart.y;}
                else {dividedVelocityIncrement.y = dividedVelocityIncrement.y - dividedVelocityStart.y;}

                //For minus values:
                if (dividedVelocityIncrement <= 0)
                {
                    if (dividedVelocityIncrement >= velocities[incrementY][incrementX]) {break;}
                }
            }
        }

        ++incrementX;
        if (incrementX == vectorWidth)
        {
            incrementX = 0;
            ++incrementY;
        }
    }
}