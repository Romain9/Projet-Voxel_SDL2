#include <iostream>
#include <chrono>
#include <SDL.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

/* Ressources utilisées

Série de vidéos sur la concéption d'un moteur 3D (utile comme initiation):
    https://www.youtube.com/watch?v=ih20l3pJoeU&list=PLrOv9FMX8xJE8NgepZR1etrsU63fDDGxO&index=22
    https://www.youtube.com/watch?v=XgMWc6LumG4&list=PLrOv9FMX8xJE8NgepZR1etrsU63fDDGxO&index=23
    https://www.youtube.com/watch?v=HXSuNxpCzdM&list=PLrOv9FMX8xJE8NgepZR1etrsU63fDDGxO&index=24
    bonus: https://www.youtube.com/watch?v=nBzCS-Y0FcY&list=PLrOv9FMX8xJE8NgepZR1etrsU63fDDGxO&index=25

Article et vidéos pour comprendre le fonctionnement de la caméra dans un moteur 3D 
(Principalement sur les 3 matrices utiles pour modéliser un objet (ici cube/voxel) et sur
les transformations nécessaire à la caméra pour visualiser cet objet):
    https://www.youtube.com/watch?v=-tonZsbHty8
    https://jsantell.com/model-view-projection/#:~:text=Model%20View%20Projection%20is%20a,which%20can%20then%20be%20rasterized.
    http://www.codinglabs.net/article_world_view_projection_matrix.aspx

Pour comprendre le bruit de Perlin (perlin-noise): 
    (L'algorithme est compliqué (en tout cas pour moi) et après toutes les recherches et l'implémentation du Perlin--Noise
    je ne maîtrise pas assez bien le sujet donc je n'ai pas osez écrire de commentaire sur les différentes fonctions
    utilisées pour éviter toutes mauvaise explications de la vidéo ci-dessous.)

    https://www.youtube.com/watch?v=22Ny_nX-YlA -> le code source entier sur le bruit de perlin viens de là

    Bonus (pour différentes manière d'approcher le perlin-noise):
        https://www.youtube.com/watch?v=6-0UaeJBumA&t=1503s
        https://www.youtube.com/watch?v=8ZEMLCnn8v0&t=698s
        https://www.youtube.com/watch?v=IKB1hWWedMk&t=1183s

Article utile globalement:
    https://blog.devgenius.io/how-i-built-a-basic-3d-graphics-engine-from-scratch-a54df82031f3

*/

struct vec2d
{
    float x = 0.0f;
    float y = 0.0f;
};




struct vec3d 
{
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    
};


struct triangle 
{
    vec3d p[3];
};

struct triangleLum
{
    triangle t;
    float lum;
};

struct mesh
{
    std::vector<triangle> tris;

    void clear()
    {
        tris.clear();
    }
};

struct mat4x4
{
    float m[4][4];
};


class engine3D
{
public:
    SDL_Event e;
    SDL_Window* window;
    SDL_Renderer* renderer;
    std::vector<SDL_FPoint> points;

public:
    engine3D(unsigned int _seed)
    {
        
        SDL_Init(SDL_INIT_VIDEO);
        window = SDL_CreateWindow("SDL Projet Voxel", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

        SDL_RenderSetScale(renderer, 2, 2);

        
        reseed(_seed);
       
    }

    void reseed(unsigned int _seed)
    {
        seed = _seed;

        // populate the permutation table in order

        for(unsigned int i = 0 ; i < 256 ; i++)
        {
            p[i] = i;
        }

        // shuffle 
        std::shuffle(std::begin(p), std::begin(p) + 256, std::default_random_engine(seed));

        // duplicate array for overflow
        for(unsigned int i = 0 ; i < 256 ; i++)
        {
            p[256 + i] = p[i];
        }
    }

    
     

    float noise2D(float x, float y)
    {
        // find smallest point of square containing target
        int xi = (int)(floorf(x)) & 255;
        int yi = (int)(floorf(y)) & 255;

        // get decimal value of each component

        x -= floorf(x);
        y -= floorf(y);

        // get smooth value from fade function (becomes weight for each dimension)

        float sx = fade(x);
        float sy = fade(y);

        // get hash value for all neighboring points

        int aa, ab, ba, bb;

        aa = p[p[xi    ] + yi    ];
        ab = p[p[xi    ] + yi + 1];
        ba = p[p[xi + 1] + yi    ];
        bb = p[p[xi + 1] + yi + 1];

        // get weighted average

        float avg = lerp(
            sy,
            lerp(   // top
                sx,
                grad(aa, x    , y    , 0),
                grad(ba, x - 1, y    , 0)
            ),
            lerp(   // bottom
                sx,
                grad(ab, x    , y - 1, 0),
                grad(bb, x - 1, y - 1, 0)
            )
        );

        return map(avg, -1, 1, 0, 1);
        
    }

    float accumulatedNoise2D(float x, float y, int octaves = 8, float lacunarity = 2.0f, float gain = 0.5f) 
    {   
        /*
            accumule un certain nombre de bruit de Perlin pour en produire un plus de variation
            (changez la valeur par défaut d'octaves pour observer un terrain avec plusieurs dunes de taille varié (octave = 8)
            ou avec 1 à 2 énormes vagues (octave = 1)).
        */
        float result = 0.0f;
        float amplitude = 1.0f;
        float frequency = 1.0f;
        float maxVal = 0.0f; // used to normalise result

        for(; octaves > 0 ; octaves--)
        {
            result += noise2D(x * frequency, y * frequency) * amplitude;

            maxVal += amplitude;

            amplitude *= gain;

            frequency *= lacunarity;
        }


        return result / maxVal;
    }

    void show(float fTime)
    {
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); // on définit la couleur de la fenêtre 
        SDL_RenderClear(renderer);  // on applique la modification appliqué au renderer (ici la couleur de fond de la fenêtre)
        
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); // le prochain ajout à la fenêtre (un point, une ligne etc...) sera blanc
    
        
        onUserUpdate(fTime);


        SDL_RenderPresent(renderer);

    }


    void drawTriangle(float x1, float y1, float x2, float y2, float x3, float y3, float lum)
    {
        SDL_RenderDrawLineF(renderer, x1, y1, x2, y2);
        // red = SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
        SDL_RenderDrawLineF(renderer, x2, y2, x3, y3);
        //  green = SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
        SDL_RenderDrawLineF(renderer, x3, y3, x1, y1);
        //  blue = SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        
        
        Uint8 r = 255;
        Uint8 g = 255;
        Uint8 b = 255;
        Uint8 a = 255;
        
        /*
            Le systeme d'éclairage est "fonctionnel" mais nécessite amélioration. 
        */
        
        getColour(lum, r, g, b, a);  // commenter pour désactiver l'éclairage


        const std::vector<SDL_Vertex> verts = // Remplis la face d'un triangle
        {
            { SDL_FPoint{ x1, y1 }, SDL_Color{ r, g, b, a }, SDL_FPoint{ 0 }, },
            { SDL_FPoint{ x2, y2 }, SDL_Color{ r, g, b, a }, SDL_FPoint{ 0 }, },
            { SDL_FPoint{ x3, y3 }, SDL_Color{ r, g, b, a }, SDL_FPoint{ 0 }, },
        };
        
        SDL_RenderGeometry(renderer, nullptr, verts.data(), verts.size(), nullptr, 0); // affiche le triangle

    }


    void getColour(float lum, Uint8 &r, Uint8 &g, Uint8 &b, Uint8 &sym)
    {

        /* 
           La couleur attribuée aux face des triangles est choisis selon
           un coefficient de luminosité (=lum).
        */

        int pixel_bw = (int)(13.0f*lum);

        int PIXEL_SOLID = 255;
        int PIXEL_THREEQUARTER = 191;
        int PIXEL_HALF = 127;
        int PIXEL_QUARTER = 63;

        switch (pixel_bw)
        {
            case 0:
            {
                r = 0;
                g = 0;
                b = 0;
                sym = PIXEL_SOLID;
                
                break;
            }
            case 1:
            {
                r = 21;
                g = 21;
                b = 21;
                sym = PIXEL_QUARTER;

                
                break;
            }
            case 2:
            {
                r = 42;
                g = 42;
                b = 42;
                sym = PIXEL_HALF;

                
                break;
            }
            case 3:
            {
                r = 63;
                g = 63;
                b = 63;
                sym = PIXEL_THREEQUARTER;

                
                break;
            }
            case 4:
            {
                r = 84;
                g = 84;
                b = 84;
                sym = PIXEL_SOLID;

                
                break;
            }
            case 5:
            {
                r = 105;
                g = 105;
                b = 105;
                sym = PIXEL_QUARTER;

                
                break;
            }
            case 6:
            {
                r = 126;
                g = 126;
                b = 126;
                sym = PIXEL_HALF;

                
                break;
            }
            case 7:
            {
                r = 147;
                g = 147;
                b = 147;
                sym = PIXEL_THREEQUARTER;

                
                break;
            }
            case 8:
            {
                r = 168;
                g = 168;
                b = 168;
                sym = PIXEL_SOLID;

                
                break;
            }
            case 9:
            {
                r = 189;
                g = 189;
                b = 189;
                sym = PIXEL_QUARTER;

                
                break;
            }
            case 10:
            {
                r = 210;
                g = 210;
                b = 210;
                sym = PIXEL_HALF;

                
                break;
            }
            case 11:
            {
                r = 231;
                g = 231;
                b = 231;
                sym = PIXEL_THREEQUARTER;

                
                break;
            }
            case 12:
            {
                r = 252;
                g = 252;
                b = 252;
                sym = PIXEL_SOLID;

                
                break;
            }
            default:
            {
                r = 0;
                g = 0;
                b = 0;
                sym = PIXEL_SOLID;

                
                break;
            }
        }
    }

    uint32_t nLehmer = 0; // Première tentative (avant le perlin-noise) pour crée une variation de hauteur entre les cubes (Un des meilleurs algorithme d'aléatoire mais le Perlin-Noise)
    uint32_t Lehmer32()   // Au final pas utilisé mais plutôt cool
    {
        nLehmer += 0xe120fc15;
        uint64_t tmp;
        tmp = (uint64_t)nLehmer * 0x4a39b70d;
        uint32_t m1 = (tmp >> 32) ^ tmp;
        tmp = (uint64_t)m1 * 0x12fad5c9;
        uint32_t m2 = (tmp >> 32) ^ tmp;
        return m2;
    }
    

    void input() // Sert uniquement pour fermer la fenêtre, j'ai préférer avoir tous les contrôles de terrain et de caméra dans la méthode OnUserUpdate.... et donc oui le code est dégueu
    {
        while(SDL_PollEvent(&e))
        {
            if(e.type == SDL_QUIT)
            {
                SDL_Quit();
                exit(0);
            }
            if(e.type == SDL_KEYDOWN)
            {
                if(e.key.keysym.sym == SDLK_ESCAPE)
                {
                    SDL_Quit();
                    exit(0);
                }
            }
        }
    }

    bool onUserCreate()
    {
        matProj = Matrix_MakeProjection(50.0f, HEIGHT / WIDTH, 0.1f, 1000.0f); 
        // creation d'une matrice de projection (matrice qui va projeté un object en 3D (ici chaque point d'un cube/voxel) sur un écran 2D)
        // Cela permet notament d'ajouter de la perspective à l'objet.

        fNoiseSeed2D = new float[nOutPutSize];
        fPerlinNoise2D = new float[nOutPutSize];

        for(int i = 0 ; i < nOutPutSize ; i++)
        {
            fNoiseSeed2D[i] = (float)rand()/(float)RAND_MAX;
        }
        
        

        return true;
    }


    bool onUserUpdate(float fElapsedTime)
    {

        matrotX = Matrix_MakeRotationX(modelRollAxis); //3.14159/3  

        matrotY = Matrix_MakeRotationY(modelYawAxis);

        matrotZ = Matrix_MakeRotationZ(modelPitchAxis);

        // Matrices qui seront utilisé pour le monde cubique

        /*
            Les matrices de rotations permettent évidement d'éffectuer une rotation sur un "objet" auquel on aurait appliqué ces matrices.
            Ici, la caméra ou le monde cubique (pour ce dérnier la variable modelPitchAxis est la seule qui est utilisé).
            Pour la caméra, les rotations possibles vont permettre une vision en première personne.
            Pour le monde cubique, la seule rotation utilisé est sur l'axe Z pour tourner le "monde" comme un plateau.
            (Note pour plus tard le monde est aussi tourné sur l'axe X pour pouvoir le visualier de profile et non de dessus,
            puisque la caméra regarde le long de l'axe Z ).
        
        */
        matTrans = Matrix_MakeTranslation(0.0f, 0.0f, 1.0f); 

        mat4x4 matWorld; 
        
        
        matWorld = Matrix_MakeIdentity();
        matWorld = Matrix_MultiplyMatrix(matWorld, matrotZ);
        matWorld = Matrix_MultiplyMatrix(matWorld, matrotX);
        matWorld = Matrix_MultiplyMatrix(matWorld, matrotY);
        matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);

        /*
            La matrice World (aussi connue comme matrice Model ou encore matrice Model To World) permet d'appliquer tous types
            de mouvement à un objet (Rotation X,Y,Z, Translation). 
            Si on veut visualiser une roue qui tourne c'est avec cette matrice que le mouvement sera appliqué.
        
        */

        vec3d vUp = { 0, 1, 0 };
        vec3d vTarget = { 0, 0, 1};
        mat4x4 matCameraRotY = Matrix_MakeRotationY(camYawAxis);
        mat4x4 matCameraRotX = Matrix_MakeRotationX(pitchAxis); 
        mat4x4 matCameraRotZ = Matrix_MakeRotationZ(rollAxis);
        mat4x4 matCameraRot = Matrix_MultiplyMatrix(matCameraRotZ, matCameraRotX);
        matCameraRot = Matrix_MultiplyMatrix(matCameraRotX, matCameraRot);
        matCameraRot = Matrix_MultiplyMatrix(matCameraRotY, matCameraRot);
        multiplyMatrixVector(vTarget, vLookDir, matCameraRot);
        vTarget = Vector_Add(vCamera, vLookDir);
        
        mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp); 

        mat4x4 matView = Matrix_QuickInverse(matCamera);

        /*
            La matrice View permet à la caméra de regarder le long d'un axe en particulier et 
            passer d'une vision de l'espace Monde dans lequel existe chaque "objet" (cube, monde cubique, caméra etc...)
            à une vision de l'espace de Vision (l'espace de caméra, grosso modo ce que voit la caméra).
        */
        
        const Uint8 *keystates = SDL_GetKeyboardState(NULL);
        
        vec3d vForward = Vector_Mul(vLookDir, 1.0f * fElapsedTime);


        /* CONTROLE CLAVIER

            xPosition et yPosition sont utilisé pour modifier le bruit de Perlin.
            Les cubes eux-mêmes ne se déplace pas, c'est une illusion d'optique en quelque sorte avec le bruit de Perlin.
            Lorsque le monde cubique est géneré, un bruit de Perlin, crée préalablement avec une seed, est appliqué aux monde cubique
            qui modifie la valeur Z (la hauteur en gros) de chaque cube.
            Ainsi lorsque on augmente ou décroit xPosition ou yPosition c'est en réalité le bruit de Perlin qu'on déplace.

            Pour une explication visuel, allez voir les lignes 670-671.

        */
        
        if(keystates[SDL_SCANCODE_D]) 
        {
            xPosition += 4.0f * fElapsedTime; 
        }
        if(keystates[SDL_SCANCODE_A])
        {
            xPosition -= 4.0f * fElapsedTime;
        }
        if(keystates[SDL_SCANCODE_W])
        {
            yPosition += 4.0f * fElapsedTime;
            
        }
        if(keystates[SDL_SCANCODE_S])
        {
            yPosition -= 4.0f * fElapsedTime;
            
        }
        if(keystates[SDL_SCANCODE_UP]) // Permet d'avancer la caméra 
        {
            
            vCamera.z += 10.0f * fElapsedTime;//Vector_Add(vCamera, vForward);
            
        }
        if(keystates[SDL_SCANCODE_DOWN]) // Permet de reculer la caméra
        {
            vCamera.z -= 10.0f * fElapsedTime; //Vector_Sub(vCamera, vForward);
            
        }
        if(keystates[SDL_SCANCODE_P]) // Fait tourner le monde cubique dans le sens des aiguilles d'une montre
        {
            modelPitchAxis -= 0.5f * fElapsedTime;
            
        }
        if(keystates[SDL_SCANCODE_O]) // Fait tourner le monde cubique dans le sens inverse des aiguilles d'une montre
        {
            modelPitchAxis += 0.5f * fElapsedTime;
            
        }
        if(keystates[SDL_SCANCODE_8]) // Bouge la caméra vers le haut
        {
            
            vCamera.y -= 4.0f;
        }
        if(keystates[SDL_SCANCODE_2]) // Bouge la caméra vers le bas
        {
            
            vCamera.y += 4.0f;
        }
        if(keystates[SDL_SCANCODE_KP_4]) // Tourne la caméra  vers la gauche (NE FONCTIONNE PAS CORRECTEMENT)
        {
            camYawAxis += 0.05f * fElapsedTime;
            
        }
        if(keystates[SDL_SCANCODE_KP_6]) // Tourne la caméra  vers la droite (NE FONCTIONNE PAS CORRECTEMENT)
        {
            
            camYawAxis -= 0.05f * fElapsedTime;
            
        }
        if(keystates[SDL_SCANCODE_KP_8]) // Tourne la caméra  vers le haut 
        {
            
            pitchAxis += sinf(0.05f * fElapsedTime);
           
        }
        if(keystates[SDL_SCANCODE_KP_2]) // Tourne la caméra  vers le bas
        {
            
            pitchAxis -= sinf(0.05f * fElapsedTime);
            
        }
        
        std::vector<triangleLum> vecTriToRaster;


        nSectorX = 60; // taille du monde cubique
        nSectorY = 60; // taille du monde cubique

        vec3d screenSector = { 0, 0, 0 };
    
        
        float fScale = 1.0f;
        float fScaleAcc = 0.0f;
        for(int x = 0 ; x < nOutPutSize ; x++)
        {
            float fNoise = 0.0f;
            for(int o = 0 ; o < 2 ; o++)
            {
                
                int nPitch = nOutPutSize >> o;
                int nSample1 = (x / nPitch) * nPitch;
                int nSample2 = (nSample1 + nPitch) % nOutPutSize;

                float fBlend = (float)(x - nSample1) / (float)nPitch;
                float fSample = (1.0f - fBlend) * fNoiseSeed2D[nSample1] + fBlend * fNoiseSeed2D[nSample2];
                fNoise += fSample * fScale;
                fScaleAcc += fScale;
                fScale = fScale / 2.0f;
            
            }
            fPerlinNoise2D[x] = fNoise / fScaleAcc;
        }

        double frequency = 8;
        double fx = 512 / frequency;
        double fy = 512 / frequency;
        
        

        auto tp1 = std::chrono::system_clock::now();

        
        
        for(screenSector.x = 0 ; screenSector.x < nSectorX ; screenSector.x++)
        {
            
            for(screenSector.y = 0 ; screenSector.y < nSectorY ; screenSector.y++)
            {

                float cptNoiseZ = (fPerlinNoise2D[(int)(screenSector.x)]*(nSectorY*2)/2.0f) + nSectorY;
                


                double pt = accumulatedNoise2D((screenSector.x+(int)xPosition) / fx, (screenSector.y+(int)yPosition) / fy); 

                //SDL_SetRenderDrawColor(renderer, pt*255, pt*255, pt*255, 255); // Lignes à décommenter pour visualiser le Perlin-Noise 
                //SDL_RenderDrawPoint(renderer, screenSector.x, screenSector.y); // Lignes à décommenter pour visualiser le Perlin-Noise 
                

                meshCube = makeCube(meshCube, screenSector.x, screenSector.y, pt*100); // cube avec la hauteur trouvée grâce au Perlin-Noise 
                

                
                for(auto tri : meshCube.tris)
                {       
                    
                
                    triTransformed.p[0] = { (tri.p[0].x - nSectorX/2)*cosf(modelPitchAxis) - (tri.p[0].y - nSectorX/2)*sinf(modelPitchAxis), (tri.p[0].x - nSectorX/2)*sinf(modelPitchAxis) + (tri.p[0].y - nSectorX/2)*cosf(modelPitchAxis), tri.p[0].z }; 
                    triTransformed.p[1] = { (tri.p[1].x - nSectorX/2)*cosf(modelPitchAxis) - (tri.p[1].y - nSectorX/2)*sinf(modelPitchAxis), (tri.p[1].x - nSectorX/2)*sinf(modelPitchAxis) + (tri.p[1].y - nSectorX/2)*cosf(modelPitchAxis), tri.p[1].z }; 
                    triTransformed.p[2] = { (tri.p[2].x - nSectorX/2)*cosf(modelPitchAxis) - (tri.p[2].y - nSectorX/2)*sinf(modelPitchAxis), (tri.p[2].x - nSectorX/2)*sinf(modelPitchAxis) + (tri.p[2].y - nSectorX/2)*cosf(modelPitchAxis), tri.p[2].z }; 
                    // rotation du terrain sur l'axe Z (formule aussi longue car on multiplie par le centre de rotation, le centre du terrain) 

                    triTransformed.p[0] = { triTransformed.p[0].x, (triTransformed.p[0].y - nSectorX/2)*cosf(3.14159/3) - (triTransformed.p[0].z - nSectorX/2)*sinf(3.14159/3), (triTransformed.p[0].y - nSectorX/2)*sinf(3.14159/3) + (triTransformed.p[0].z - nSectorX/2)*cosf(3.14159/3)}; 
                    triTransformed.p[1] = { triTransformed.p[1].x, (triTransformed.p[1].y - nSectorX/2)*cosf(3.14159/3) - (triTransformed.p[1].z - nSectorX/2)*sinf(3.14159/3), (triTransformed.p[1].y - nSectorX/2)*sinf(3.14159/3) + (triTransformed.p[1].z - nSectorX/2)*cosf(3.14159/3)}; 
                    triTransformed.p[2] = { triTransformed.p[2].x, (triTransformed.p[2].y - nSectorX/2)*cosf(3.14159/3) - (triTransformed.p[2].z - nSectorX/2)*sinf(3.14159/3), (triTransformed.p[2].y - nSectorX/2)*sinf(3.14159/3) + (triTransformed.p[2].z - nSectorX/2)*cosf(3.14159/3)}; 
                    // rotation du terrain sur l'axe X (même réflexion pour la longueur de la formule)


                    
                    line1 = Vector_Sub(triTransformed.p[1], triTransformed.p[0]);
                    line2 = Vector_Sub(triTransformed.p[2], triTransformed.p[0]);

                    normal = Vector_CrossProduct(line1, line2); // cross product = produit véctoriel

                    normal = Vector_Normalise(normal); // permet de trouver la diréction d'un vecteur

                    /*

                        Les normals sont des vecteurs qui sont perpendiculaire à une surface, cela va permettre, grossièrement, à notre code
                        de savoir quel face d'un cube est cachée par une autre en vérifiant leurs direction. 

                    */


                    vec3d vCameraRay = Vector_Sub(triTransformed.p[0], vCamera);



                    if(Vector_DotProduct(normal, vCameraRay) < 0.0f) // On affiche uniquement les faces qui oppose la caméra, donc visible à la caméra.
                    {   
                    

                        lightDirection = Vector_Normalise(lightDirection); // Dans quel direction viens la lumière

                        float lightDotProduct = std::max(0.1f, Vector_DotProduct(lightDirection, normal));

                        multiplyMatrixVector(triTransformed.p[0], triViewed.p[0], matView); // on passe les triangles du cube de l'espace model à l'espace de vision (ce que voit la caméra)
                        multiplyMatrixVector(triTransformed.p[1], triViewed.p[1], matView); // on passe les triangles du cube de l'espace model à l'espace de vision (ce que voit la caméra)
                        multiplyMatrixVector(triTransformed.p[2], triViewed.p[2], matView); // on passe les triangles du cube de l'espace model à l'espace de vision (ce que voit la caméra)

                        int nClippedTriangles = 0;
                        triangle clipped[2];
                        nClippedTriangles = Triangle_ClipAgainstPlane({0.0f, 0.0f, 0.1f}, {0.0f, 0.0f, 1.0f}, triViewed, clipped[0], clipped[1]); 
                        // Les trois lignes plus hautes sont une tentatives de clipping (en gros parmis les triangles qui sont en dehors du champ de vision ne seront pas affichés)
                        //...ça marche pas

                        
                        
                        for(int n = 0 ; n < nClippedTriangles ; n++) // on affiche les triangles dans le champs de vision
                        {
                           
                            multiplyMatrixVector(clipped[n].p[0], triProj.p[0], matProj); // on applique la matrice de projection (ce qui permet d'ajouté la perspective)
                            multiplyMatrixVector(clipped[n].p[1], triProj.p[1], matProj); // on applique la matrice de projection (ce qui permet d'ajouté la perspective)
                            multiplyMatrixVector(clipped[n].p[2], triProj.p[2], matProj); // on applique la matrice de projection (ce qui permet d'ajouté la perspective)


                            triProj.p[0].x *= 0.5f * 640.0f;
                            triProj.p[1].x *= 0.5f * 640.0f;
                            triProj.p[2].x *= 0.5f * 640.0f;
                            triProj.p[0].y *= 0.5f * 360.0f;
                            triProj.p[1].y *= 0.5f * 360.0f;
                            triProj.p[2].y *= 0.5f * 360.0f;
                            // on éloigne le cube de la caméra

                            /* Voir prochain commentaire
                                triangleLum tProj = {triProj, lightDotProduct};
                                vecTriToRaster.push_back(tProj);
                                vecTriToRaster[0].lum = lightDotProduct;
                            */


                            drawTriangle(triProj.p[0].x, triProj.p[0].y,
                            triProj.p[1].x, triProj.p[1].y,
                            triProj.p[2].x, triProj.p[2].y, lightDotProduct);
                            
                            
                        } 
                        
                    }
                }

                /* Pas utilisé pour des raison de performance (en gros ça permet de trier selon quel face est plus proche et donc d'afficher en premier ces faces)

                sort(vecTriToRaster.begin(), vecTriToRaster.end(), [](triangleLum t1, triangleLum t2)
                {
                    float z1 = (t1.t.p[0].z + t1.t.p[1].z + t1.t.p[2].z) / 3.0f;
                    float z2 = (t2.t.p[0].z + t2.t.p[1].z + t2.t.p[2].z) / 3.0f;
                    return z1 > z2;
                });
                
                for(auto &triProjected : vecTriToRaster)
                {
                    drawTriangle(triProjected.t.p[0].x, triProjected.t.p[0].y,
                            triProjected.t.p[1].x, triProjected.t.p[1].y,
                            triProjected.t.p[2].x, triProjected.t.p[2].y, triProjected.lum);
                }
                */
                
            }
            
            
            
        }
        
       
        
        
        
        auto tp2 = std::chrono::system_clock::now();

        elapsedTime = tp2 - tp1; // pratique pour débug, calcule le temps entre tp1 et tp2...
        
        if(keystates[SDL_SCANCODE_SPACE])
        {
            std::cout << std::to_string(elapsedTime.count()) << std::endl; //...l'affiche d'une pression de la barre espace
        }
        

        return true;
    }

private:
    const float WIDTH = 1280.0f;
    const float HEIGHT = 720.0f;

    int nSectorX = 8;
    int nSectorY = 8;

    float floatArr[5] = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f };

    int nOutPutSize = 80*80;

    float *fNoiseSeed2D = nullptr;
    float *fPerlinNoise2D = nullptr;

    std::chrono::duration<float> elapsedTime;

    mesh meshCube;

    mat4x4 matrotX, matrotY, matrotZ, matTrans;
    mat4x4 matProj;

    mat4x4 modToWorld;
    mat4x4 worldToView;
    mat4x4 matOrthProj;
    vec3d normal, line1, line2;

    triangle triProj, triOrthProj, triTranslated, triViewed,  triRotatedZ, triRotatedX, triRotatedY, triTransformed;
    
    vec3d vCamera = { -100.0f, -200.0f, -200.0f };
    vec3d vLookDir;

    float xPosition = 0.0f;
    float yPosition = 0.0f;

    float camYawAxis = 0.0f;
    float pitchAxis = 0.0f;
    float rollAxis = 0.0f;

    float modelYawAxis = 0.0f;
    float modelRollAxis = 0.0f;
    float modelPitchAxis = 0.0f;

    vec3d lightDirection = { 0.0f, 0.0f, -1.0f };

    float lightDotProduct;

    float fTheta = 0.0f;

    // perlin noise variables

    unsigned char p[512];   // permutation table

    
    unsigned int seed; // current seed



    // fade function for perlin noise
    // https://www.wolframalpha.com/input/?i=plot+6*t%5E5+-+15*t%5E4+%2B+10*t%5E3+from+0+to+1   Pour visualiser la fonction

    inline float fade(float t)
    {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }

    // interpolation linéaire (De ce que j'ai compris, cela permet d'estimer où une donnée pourrait existé en dehors des valeurs connue)

    inline float lerp(float t, float a, float b)
    {
        return a + t * (b - a);
    }

    // map function

    inline float map(float val, float ogMin, float ogMax, float newMin, float newMax)
    {
        // get original range proportion

        float prop = (val - ogMin) / (ogMax - ogMin);

        return lerp(prop, newMin, newMax);
    }

    // gradient function -> calculate dot product between gradient vector and distance vector

    inline float grad(int hash, float x, float y, float z)
    {
        // convert last 4 bit of the hash into one of 12 possible gradients

        int h = hash & 0b1111; // = hash & 15 = hash % 16 = get last 4 bits

        // if 1st bit is 1, set to x, otherwise set to y

        float u = h < 0b1000 ? x : y;

        // if 1st/2nd bits are 0, set to y
        // if 1st/2nd bits are 1, set to x
        // else set to z

        float v = h < 0b0100 ? y : h == 0b1100 || h == 0b1110 ? x : z;

        // make u, v negative based on the last 2 bits, then add results
        // adding results is like the dot product because the gradient components are 1,
        // so the results of the dot product is adding the distance components

        return ((h & 0b0001) == 0 ? u : -u) + ((h & 0b0010) == 0 ? v : -v);
    }

    mesh makeCube(mesh cube, float x, float y, float z)
    {
        cube.tris = {
        
            // SOUTH
            { x, y, z,    x, y+1, z,    x+1, y+1, z },
            { x, y, z,    x+1, y+1, z,    x+1, y, z },

            // EAST                                                      
            { x+1, y, z,    x+1, y+1, z,    x+1, y+1, z+1 },
            { x+1, y, z,    x+1, y+1, z+1,    x+1, y, z+1 },

            // NORTH                                                     
            { x+1, y, z+1,    x+1, y+1, z+1,    x, y+1, z+1 },
            { x+1, y, z+1,    x, y+1, z+1,    x, y, z+1 },

            // WEST                                                      
            { x, y, z+1,    x, y+1, z+1,    x, y+1, z },
            { x, y, z+1,    x, y+1, z,    x, y, z },

            // TOP                                                       
            { x, y+1, z,    x, y+1, z+1,    x+1, y+1, z+1 },
            { x, y+1, z,    x+1, y+1, z+1,    x+1, y+1, z },

		    // BOTTOM                                                    
	        { x+1, y, z+1,    x, y, z+1,    x, y, z },
            { x+1, y, z+1,    x, y, z,    x+1, y, z },
        };

       return cube; 
    }

    void multiplyMatrixVector(vec3d &i, vec3d &o, mat4x4 &m)
    {
        o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
		o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
		o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
		float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];


		if (w != 0.0f)
		{
			o.x /= w; o.y /= w; o.z /= w;
		}
        
    }


    mat4x4 Matrix_MakeIdentity() // Matrice identité, algèbre linéaire twa meme tu C
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[0][1] = 0.0f;
        matrix.m[0][2] = 0.0f;
        matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = 0.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[1][2] = 0.0f;
        matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = 0.0f;
        matrix.m[2][1] = 0.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = 0.0f;
        matrix.m[3][1] = 0.0f;
        matrix.m[3][2] = 0.0f;
        matrix.m[3][3] = 1.0f;

        return matrix;
    }

    mat4x4 Matrix_MakeRotationX(float fAngleRad) // La matrice de base pour toute rotation sur l'axe X (check wikipedia si tu me crois pas)
    {
       mat4x4 matrix;
       matrix.m[0][0] = 1.0f;
       matrix.m[0][1] = 0.0f;
       matrix.m[0][2] = 0.0f;
       matrix.m[0][3] = 0.0f;
       matrix.m[1][0] = 0.0f;
       matrix.m[1][1] = cosf(fAngleRad);
       matrix.m[1][2] = sinf(fAngleRad);
       matrix.m[1][3] = 0.0f;
       matrix.m[2][0] = 0.0f;
       matrix.m[2][1] = -sinf(fAngleRad);
       matrix.m[2][2] = cosf(fAngleRad);
       matrix.m[2][3] = 0.0f;
       matrix.m[3][0] = 0.0f;
       matrix.m[3][1] = 0.0f;
       matrix.m[3][2] = 0.0f;
       matrix.m[3][3] = 1.0f;

       return matrix;
    }

    mat4x4 Matrix_MakeRotationY(float fAngleRad) // La matrice de base pour toute rotation sur l'axe Y
    {
       mat4x4 matrix;
       matrix.m[0][0] = cosf(fAngleRad);
       matrix.m[0][1] = 0.0f;
       matrix.m[0][2] = sinf(fAngleRad);
       matrix.m[0][3] = 0.0f;
       matrix.m[1][0] = 0.0f;
       matrix.m[1][1] = 1.0f;
       matrix.m[1][2] = 0.0f;
       matrix.m[1][3] = 0.0f;
       matrix.m[2][0] = -sinf(fAngleRad);
       matrix.m[2][1] = 0.0f;
       matrix.m[2][2] = cosf(fAngleRad);
       matrix.m[2][3] = 0.0f;
       matrix.m[3][0] = 0.0f;
       matrix.m[3][1] = 0.0f;
       matrix.m[3][2] = 0.0f;
       matrix.m[3][3] = 1.0f;
       
       return matrix;
    }


    mat4x4 Matrix_MakeRotationZ(float fAngleRad) // La matrice de base pour toute rotation sur l'axe Z
    {
       mat4x4 matrix;
       matrix.m[0][0] = cosf(fAngleRad);
       matrix.m[0][1] = sinf(fAngleRad);
       matrix.m[0][2] = 0.0f;
       matrix.m[0][3] = 0.0f;
       matrix.m[1][0] = -sinf(fAngleRad);
       matrix.m[1][1] = cosf(fAngleRad);
       matrix.m[1][2] = 0.0f;
       matrix.m[1][3] = 0.0f;
       matrix.m[2][0] = 0.0f;
       matrix.m[2][1] = 0.0f;
       matrix.m[2][2] = 1.0f;
       matrix.m[2][3] = 0.0f;
       matrix.m[3][0] = 0.0f;
       matrix.m[3][1] = 0.0f;
       matrix.m[3][2] = 0.0f;
       matrix.m[3][3] = 0.0f;
       
       return matrix;
    }

    

    mat4x4 Matrix_MakeTranslation(float x, float y, float z) // Matrice de translation, prise sur la 3ème vidéo des ressources (tout en haut)
    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[0][1] = 0.0f;
        matrix.m[0][2] = 0.0f;
        matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = 0.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[1][2] = 0.0f;
        matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = 0.0f;
        matrix.m[2][1] = 0.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = x;
        matrix.m[3][1] = y;
        matrix.m[3][2] = z;
        matrix.m[3][3] = 1.0f;

        return matrix;
    }

    mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar) // Matrice de projection, tout comme le dernier commentaire
    {
        float fFovRad = 1.0f / tanf((fFovDegrees * 3.14159) / (180.0f * 2)); // FOV = Field of view; ici en radian
        mat4x4 matrix;
        matrix.m[0][0] = fAspectRatio * fFovRad;
        matrix.m[0][1] = 0.0f;
        matrix.m[0][2] = 0.0f;
        matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = 0.0f;
        matrix.m[1][1] = fFovRad;
        matrix.m[1][2] = 0.0f;
        matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = 0.0f;
        matrix.m[2][1] = 0.0f;
        matrix.m[2][2] = fFar / (fFar - fNear);
        matrix.m[2][3] = 1.0f;
        matrix.m[3][0] = 0.0f;
        matrix.m[3][1] = 0.0f;
        matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
        matrix.m[3][3] = 0.0f;

        return matrix;
    }

    
    mat4x4 Matrix_MultiplyMatrix(mat4x4 &m1, mat4x4 &m2) // Formule pour multiplier une matrice avec une autre
    {
        mat4x4 matrix;
        for(int c = 0; c < 4 ; c++)
        {
            for(int r = 0 ; r < 4 ; r++)
            {
                matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
            }
        }

        return matrix;
    }

    mat4x4 Matrix_PointAt(vec3d &pos, vec3d &target, vec3d &up) // Matrice utile pour la caméra qui permet de choisir sur quel axe la caméra regarde
    {
        vec3d newForward = Vector_Sub(target, pos);
        newForward = Vector_Normalise(newForward);

        vec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
        vec3d newUp = Vector_Sub(up, a);
        newUp = Vector_Normalise(newUp);

        vec3d newRight = Vector_CrossProduct(newUp, newForward);

        mat4x4 matrix;
        matrix.m[0][0] = newRight.x;
        matrix.m[0][1] = newRight.y;
        matrix.m[0][2] = newRight.z;
        matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = newUp.x;
        matrix.m[1][1] = newUp.y;
        matrix.m[1][2] = newUp.z;
        matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = newForward.x;
        matrix.m[2][1] = newForward.y;
        matrix.m[2][2] = newForward.z;
        matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = pos.x,
        matrix.m[3][1] = pos.y;
        matrix.m[3][2] = pos.z;
        matrix.m[3][3] = 1.0f;
        return matrix;
    }

    mat4x4 Matrix_QuickInverse(mat4x4 &m) // Fonction pour inverser une matrice (cette fonction est spécifique à la matrice PointAt)
    {
        mat4x4 matrix;
        matrix.m[0][0] =  m.m[0][0];
        matrix.m[0][1] =  m.m[1][0];
        matrix.m[0][2] =  m.m[2][0];
        matrix.m[0][3] =  0.0f;
        matrix.m[1][0] =  m.m[0][1];
        matrix.m[1][1] =  m.m[1][1];
        matrix.m[1][2] =  m.m[2][1];
        matrix.m[1][3] =  0.0f;
        matrix.m[2][0] =  m.m[0][2];
        matrix.m[2][1] =  m.m[1][2];
        matrix.m[2][2] =  m.m[2][2];
        matrix.m[2][3] =  0.0f;
        matrix.m[3][0] =  -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
        matrix.m[3][1] =  -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
        matrix.m[3][2] =  -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
        matrix.m[3][3] =  1.0f;
        return matrix;

    }

    vec3d Vector_Add(vec3d &v1, vec3d &v2)
    {
        return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    }

    vec3d Vector_Sub(vec3d &v1, vec3d &v2)
    {
        return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    }

    vec3d Vector_Mul(vec3d &v1, float k)
    {
        return { v1.x * k, v1.y * k, v1.z * k };
    }

    vec3d Vector_Div(vec3d &v1, float k)
    {
        return { v1.x / k, v1.y / k, v1.z / k };
    }

    float Vector_DotProduct(vec3d &v1, vec3d &v2)
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    float Vector_Length(vec3d &v)
    {
        return sqrtf(Vector_DotProduct(v, v));
    }

    vec3d Vector_Normalise(vec3d &v)
    {
        float l = Vector_Length(v);
        return { v.x / l, v.y / l, v.z / l };
    }

    vec3d Vector_CrossProduct(vec3d &v1, vec3d &v2)
    {
        vec3d v;
        v.x = v1.y * v2.z - v1.z * v2.y;
        v.y = v1.z * v2.x - v1.x * v2.z;
        v.z = v1.x * v2.y - v1.y * v2.x;
        return v;
    }

    vec3d Vector_IntersectPlane(vec3d &plane_p, vec3d &plane_n, vec3d &lineStart, vec3d &lineEnd) // Permet de savoir quel triangle entre en contact avec le plan de vision
    {
       plane_n = Vector_Normalise(plane_n);
       float plane_d = -Vector_DotProduct(plane_n, plane_p);
       float ad = Vector_DotProduct(lineStart, plane_n);
       float bd = Vector_DotProduct(lineEnd, plane_n);
       float t = (-plane_d - ad) / (bd - ad);
       vec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
       vec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
       return Vector_Add(lineStart, lineToIntersect);

    }

    int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle &in_tri, triangle &out_tri1, triangle &out_tri2) // Permet de savoir quel triangle va être clippé par le plan de vision
    {
        plane_n = Vector_Normalise(plane_n);

        auto dist = [&](vec3d &p)
        {
            vec3d n = Vector_Normalise(p);
            return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
        };

        vec3d* inside_points[3];
        int nInsidePointCount = 0;
        
        vec3d* outside_points[3];
        int nOutsidePointCount = 0;

        float d0 = dist(in_tri.p[0]);
        float d1 = dist(in_tri.p[1]);
        float d2 = dist(in_tri.p[2]);

        if(d0 >= 0)
        {
            inside_points[nInsidePointCount++] = &in_tri.p[0];
        }
        else
        {
            outside_points[nOutsidePointCount++] = &in_tri.p[0];
        }

        if(d1 >= 0)
        {
            inside_points[nInsidePointCount++] = &in_tri.p[1];
        }
        else
        {
            outside_points[nOutsidePointCount++] = &in_tri.p[1];
        }

        if(d2 >= 0)
        {
            inside_points[nInsidePointCount++] = &in_tri.p[2];
        }
        else
        {
            outside_points[nOutsidePointCount++] = &in_tri.p[2];
        }


        if(nInsidePointCount == 0)
        {
            return 0;
        }

        if(nInsidePointCount == 3)
        {
            out_tri1 = in_tri;

            return 1;
        }

        if(nInsidePointCount == 1 && nOutsidePointCount == 2)
        {
            out_tri1.p[0] = *inside_points[0];

            out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
            out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

            return 1;
        }

        if(nInsidePointCount == 2 && nOutsidePointCount == 1)
        {
            out_tri1.p[0] = *inside_points[0];

            out_tri1.p[1] = *inside_points[1];

            out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

            out_tri2.p[0] = *inside_points[1];

            out_tri2.p[1] = out_tri1.p[2];

            out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

            return 2;
        }

        return 0;
    }

};


int main(int argc, char* argv[])
{
    engine3D demo(1564865545); // Seed pour le perlin-noise, pour identifier un "monde" (dans le contexte du projet)
    
    demo.onUserCreate(); // crée une fenêtre avec SDL2
    float currentTime = SDL_GetTicks() / 1000.0f;  
    /*
        Le "temps", utilisé principalement pour calculer la vitesse au mouvement de caméra (cf:controle clavier après la ligne 532)
        En physique Vitesse = distance / temps. Ici les distances sont fixe et on multiplie par le temps car on divise le temps par 1000. 
    */
    
    while(true)
    {
        
        float startTime;
        demo.show(currentTime);
        demo.input();
        if(currentTime - startTime > 50.0f) // On vérifie que le "temps" n'est pas trop élevé. 
        {
            startTime = SDL_GetTicks() / 1000.0f; // On "reset" plus ou moins le "temps" pour éviter que les mouvements de caméra 
                                                  // deviennent trop rapide.
        }
    }
    

    return 0;
}
// "La programmation modulaire ? Connais pas..." - R0m1 2022-2023