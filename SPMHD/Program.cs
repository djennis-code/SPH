using Raylib_cs;
using System.Numerics;

namespace SPH;

internal static class Program
{
    public static Vector2 domainSize = new(1.5f, 1.0f); //  meter

    public static float upScale = 720f; //  pixels / meter
    public static Vector2 windowSize = upScale * domainSize;

    public static float currTime = 0f;

    public static HydroDynamics simulation = new();

    public static void InitSPH()
    {
        simulation.particleCount = 4900; // perfect squares are nice to work with
        simulation.smoothingRadius = 4e-2f; //  meter
        simulation.NormaliseKernels();

        simulation.domainThickness = 1f;
        simulation.wallForce = 1e5f; //  meter / second^2; rigidity of the walls
        simulation.gridSpacing = 0.2f * simulation.smoothingRadius; //  meter
        simulation.gravity = 9.81f; //  meter / second^2

        simulation.gamma = 7f; //  exponent
        simulation.refDensity = 1e3f; //  Kg / meter^3
        simulation.cSound = 1.5e2f; //  meter / second
        simulation.pressureConst = simulation.refDensity * simulation.cSound * simulation.cSound / simulation.gamma; // Pascal; B = rho * cs^2 / gamma
        simulation.pressureExt = 1e5f; //  Pascal

        simulation.particleMass = simulation.refDensity * simulation.gridSpacing * simulation.gridSpacing; //  Kg; m = rho * V
        simulation.particleRadius = 4.0f; //  pixels    
        simulation.bounds = domainSize;
        simulation.spawnOffset = new Vector2(0.0f, 0.0f);

        simulation.deltaTime = 1f * simulation.smoothingRadius / simulation.cSound; //  seconds / frame
        simulation.viscosity = 1e-1f; //  Pascal * second

        simulation.CreateBuffers();
        simulation.InitParticles();
        simulation.InitHashMap();
    }

    public static void Main()
    {
        InitSPH();

        Raylib.InitWindow((int)windowSize.X, (int)windowSize.Y, "Smoothed Particle Magneto Dynamics");

        while (!Raylib.WindowShouldClose())
        {
            Raylib.BeginDrawing();

            Raylib.ClearBackground(Color.Black);

            simulation.IntegrateEuler();
            
            RenderDensity(5);
            //RenderParticles();

            Raylib.DrawFPS(5, 5);
            Raylib.EndDrawing();

            currTime += simulation.deltaTime;
            Console.Write("Time (s): ");
            Console.WriteLine(currTime);
        }

        Raylib.CloseWindow();
    }

    public static void RenderParticles()
    {
        int count = simulation.particleCount;
        float radius = simulation.particleRadius;
        Vector2[] positions = simulation.positions;

        for (uint i = 0; i < count; i++)
        {
            Color col = Color.White;
            Vector2 pos = upScale * positions[i];

            Raylib.DrawPixel((int)pos.X, (int)pos.Y, col);
        }
    }

    public static void RenderDensity(int granularity)
    {
        float fac = 1.0f / upScale;
        float refDensity = simulation.refDensity;

        int dx = granularity;
        int dy = granularity;
        int screenWidth = (int) windowSize.X / dx;
        int screenHeight = (int) windowSize.Y / dy;
        int n = screenWidth * screenHeight;

        float[] densityBuffer = new float[n];

        Parallel.For(0, screenWidth, x =>
        {
            Parallel.For(0, screenHeight, y =>
            {
                int index = y * screenWidth + x;
                Vector2 samplePosition = new Vector2(x * dx, y * dy) * fac;

                densityBuffer[index] = simulation.ComputeDensityAtPoint(samplePosition);
            });
        });

        for (int i = 0; i < n; i++)
        {
            float density = densityBuffer[i];

            int y = (int) (i / screenWidth) * dy;
            int x = (i % screenWidth) * dx;


            Color col = DensityToColor(density, refDensity);
            Raylib.DrawRectangle(x, y, dx, dy, col);
        }

        // Print total density to see whether continuity relation holds true
        //float totalDensity = densityBuffer.Sum();
        //Console.Write("Total Density: ");
        //Console.WriteLine(totalDensity);
        // Answer: yes (close enough, may fluctuate by 0.1%)
    }

    public static Color DensityToColor(float density, float refDensity)
    {
        float fac = density / refDensity;
        int b = (int) (fac * 127);
        return new Color(0, 0, b, 255);
    }
}

