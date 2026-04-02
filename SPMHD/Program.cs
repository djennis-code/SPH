using Raylib_cs;
using System.Numerics;
using System.Runtime.InteropServices;

namespace SPH;

internal static class Program
{
    public static Vector2 domainSize = new(1.0f, 1.0f); //  meter

    public static float upScale = 800f; //  pixels / meter
    public static Vector2 windowSize = upScale * domainSize;

    public static float currTime = 0f;

    public static HydroDynamics simulation = new();

    public static void InitSPH()
    {
        simulation.particleCount = 2500; // perfect squares are nice to work with
        simulation.smoothingRadius = 4e-2f; //  meter
        simulation.NormaliseKernels();

        simulation.domainThickness = 1f;
        simulation.wallForce = 1e5f; //  meter / second^2; rigidity of the walls
        simulation.gridSpacing = 0.2f * simulation.smoothingRadius; //  meter
        simulation.gravity = 0f; //  meter / second^2

        simulation.gamma = 7f; //  exponent
        simulation.refDensity = 1e3f; //  Kg / meter^3
        simulation.cSound = 1.5e2f; //  meter / second
        simulation.pressureConst = simulation.refDensity * simulation.cSound * simulation.cSound / simulation.gamma; // Pascal; K = rho * cs^2 / gamma
        simulation.pressureExt = 1e5f; //  Pascal

        simulation.relMu = 1e0f; // relative permeability
        simulation.Bstrength = 0f; // Tesla
        simulation.Bdir = new(0f, 1f); // direction of the magnetic field

        simulation.particleMass = simulation.refDensity * simulation.gridSpacing * simulation.gridSpacing; //  Kg; m = rho * V 
        simulation.bounds = domainSize;
        simulation.spawnOffset = new(0.0f, 0.0f); // position of the fluids center of mass

        simulation.deltaTime = 0.2f * simulation.smoothingRadius / simulation.cSound; //  seconds / frame
        simulation.viscosity = 1e0f; //  Pascal * second

        simulation.CreateBuffers();
        simulation.InitParticles();
        simulation.InitHashMap();
    }

    public static void Main()
    {
        InitSPH();

        int gridSize = (int) (upScale * simulation.smoothingRadius / 2);

        Raylib.InitWindow((int)windowSize.X, (int)windowSize.Y, "Smoothed Particle Magneto Dynamics");

        simulation.IntegrateVerlet();

        while (!Raylib.WindowShouldClose())
        {
            Raylib.BeginDrawing();

            Raylib.ClearBackground(Color.Black);

            simulation.IntegrateVerlet();

            //RenderDensity(gridSize);
            RenderDivergenceB(gridSize);
            RenderParticles();

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

                densityBuffer[index] = simulation.SampleDensity(samplePosition);
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
    }

    public static Color DensityToColor(float density, float refDensity)
    {
        Color low = new(0, 0, 255, 255);
        Color mid = new(0, 0, 0, 255);
        Color high = new(255, 0, 0, 255);

        float fac = (density / refDensity) - 1f;

        if (fac < 0) return Raylib.ColorLerp(mid, low, -fac);
        else return Raylib.ColorLerp(mid, high, fac);
    }

    public static void RenderDivergenceB(int granularity)
    {
        float fac = 1.0f / upScale;

        int dx = granularity;
        int dy = granularity;
        int screenWidth = (int)windowSize.X / dx;
        int screenHeight = (int)windowSize.Y / dy;
        int n = screenWidth * screenHeight;

        float[] divergenceBuffer = new float[n];

        Parallel.For(0, screenWidth, x =>
        {
            Parallel.For(0, screenHeight, y =>
            {
                int index = y * screenWidth + x;
                Vector2 samplePosition = new Vector2(x * dx, y * dy) * fac;

                divergenceBuffer[index] = simulation.SampleDivB(samplePosition);
            });
        });

        for (int i = 0; i < n; i++)
        {
            float div = divergenceBuffer[i];

            int y = (int)(i / screenWidth) * dy;
            int x = (i % screenWidth) * dx;


            Color col = DivergenceToColor(div);
            Raylib.DrawRectangle(x, y, dx, dy, col);
        }
    }

    public static Color DivergenceToColor(float div)
    {
        Color low = new(0, 0, 255, 255);
        Color mid = new(0, 0, 0, 255);
        Color high = new(255, 0, 0, 255);

        float B = simulation.relMu / simulation.Bstrength;
        float fac = Math.Clamp(1e-2f * div / B, -1f, 1f);

        if (fac < 0) return Raylib.ColorLerp(mid, low, -fac);
        else return Raylib.ColorLerp(mid, high, fac);
    }
}

