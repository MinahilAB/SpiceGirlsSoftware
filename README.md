# SpiceGirlsSoftware
#serial diagram
```mermaid
graph TD
    %% Main flow
    A[Start: atmosphere_model] --> B[Initialize System]
    B --> C[Call init subroutine]

    %% Initialization
    C --> C1[Allocate Data Structures]
    C1 --> C2[Initialize variables: oldstat, newstat, flux, tend, ref]
    C2 --> C3[Compute grid spacing and timestep: dx, dz, dt]
    C3 --> C4[Initialize Thermal State]
    C4 --> C5[Set Reference State]

    C5 --> D[Compute total_mass_energy: mass0, te0]
    D --> E[Create output files]
    E --> F[Write initial state record]

    %% Time loop
    F --> G{Time Loop: etime < sim_time?}
    
    G -->|Yes| H[Runge-Kutta Integration]

    %% Runge-Kutta directional steps
    H --> H1{dimswitch?}
    H1 -->|True| H2[X-direction RK steps]
    H1 -->|False| H3[Z-direction RK steps]

    H2 --> H2a[step 1: dt/3, DIR_X]
    H2a --> H2b[step 2: dt/2, DIR_X]
    H2b --> H2c[step 3: dt, DIR_X]
    H2c --> H4[Z-direction RK steps]

    H3 --> H3a[step 1: dt/3, DIR_Z]
    H3a --> H3b[step 2: dt/2, DIR_Z]
    H3b --> H3c[step 3: dt, DIR_Z]
    H3c --> H5[X-direction RK steps]

    H4 --> H4a[step 1: dt/3, DIR_Z]
    H4a --> H4b[step 2: dt/2, DIR_Z]
    H4b --> H4c[step 3: dt, DIR_Z]
    H4c --> H6[Toggle dimswitch]

    H5 --> H5a[step 1: dt/3, DIR_X]
    H5a --> H5b[step 2: dt/2, DIR_X]
    H5b --> H5c[step 3: dt, DIR_X]
    H5c --> H6

    H6 --> I[Update simulation time and output counter]

    %% Output
    I --> J{Output counter >= output frequency?}
    J -->|Yes| K[Write output record]
    J -->|No| G
    K --> G

    %% End of simulation
    G -->|No| L[Compute total_mass_energy: mass1, te1]
    L --> M[Close output files]
    M --> N[Print conservation diagnostics]
    N --> O[Finalize]
    O --> P[End Simulation]

    %% Soft pastel styling with black text
    style A fill:#d0ebff,stroke:#000,stroke-width:1px,color:#000
    style P fill:#d0ebff,stroke:#000,stroke-width:1px,color:#000
    style C fill:#e6d0ff,stroke:#9b59b6,color:#000
    style H fill:#fff1d6,stroke:#d4a017,color:#000
    style G fill:#ffe5e5,stroke:#e74c3c,color:#000
    style J fill:#ffe5e5,stroke:#e74c3c,color:#000
    style D fill:#d0f0d6,stroke:#27ae60,color:#000
    style F fill:#fff8d6,stroke:#f1c40f,color:#000
    style M fill:#d0f0d6,stroke:#27ae60,color:#000
    style N fill:#fff8d6,stroke:#f39c12,color:#000

    %% Step Details subgraph
    subgraph "Step Subroutine Details"
        S1[step subroutine] --> S2{Direction?}
        S2 -->|DIR_X| S3[xtend]
        S2 -->|DIR_Z| S4[ztend]
        S3 --> S5[exchange_halo_x]
        S4 --> S6[exchange_halo_z]
        S5 --> S7[Compute Flux]
        S6 --> S7
        S7 --> S8[Compute Tendencies]
        S8 --> S9[Update state: s2 = s0 + dt*tend]
        style S1 fill:#f0f0f0,color:#000
        style S2 fill:#f0f0f0,color:#000
        style S3 fill:#fdf1e5,color:#000
        style S4 fill:#fdf1e5,color:#000
        style S5 fill:#fdf1e5,color:#000
        style S6 fill:#fdf1e5,color:#000
        style S7 fill:#fff8d6,color:#000
        style S8 fill:#fff8d6,color:#000
        style S9 fill:#fff8d6,color:#000
    end

    %% Data types subgraph
    subgraph "Data Types"
        D1[atmospheric_state] --> C1
        D2[atmospheric_flux] --> S7
        D3[atmospheric_tendency] --> S8
        D4[reference_state] --> C5
        style D1 fill:#e0f7fa,color:#000
        style D2 fill:#e0f7fa,color:#000
        style D3 fill:#e0f7fa,color:#000
        style D4 fill:#e0f7fa,color:#000
    end
```
# SpiceGirlsSoftware
#parallel diagram
```mermaid

  flowchart TD
    Start([Program Start]) --> InitMPI[MPI_Init: Initialize MPI]
    InitMPI --> GetRank[Get MPI rank & size]
    GetRank --> CreateCart[MPI_Cart_create: Create 1D Cartesian topology]
    CreateCart --> GetNeighbors[MPI_Cart_shift: Get left/right neighbor ranks]
    
    GetNeighbors --> CheckGPU{OpenACC enabled?}
    CheckGPU -->|Yes| DetectGPU[acc_get_num_devices: Detect NVIDIA GPUs]
    DetectGPU --> AssignGPU[acc_set_device_num: Assign GPU to rank]
    CheckGPU -->|No| DomainDecomp
    
    AssignGPU --> DomainDecomp[Domain Decomposition:<br/>Divide nx by MPI ranks]
    DomainDecomp --> AllocMem[Allocate memory for:<br/>oldstat, newstat, flux, tend, ref]
    
    AllocMem --> AccData{OpenACC enabled?}
    AccData -->|Yes| CopyToGPU[!$acc enter data:<br/>Copy state arrays to GPU]
    AccData -->|No| AllocBuffers
    CopyToGPU --> AllocBuffers
    
    AllocBuffers --> AllocBuf{OpenACC enabled?}
    AllocBuf -->|Yes| AllocGPUBuf[Allocate GPU buffers:<br/>send_left_d, recv_left_d, etc.]
    AllocBuf -->|No| AllocCPUBuf[Allocate CPU buffers:<br/>send_left, recv_left, etc.]
    
    AllocGPUBuf --> MainLoop
    AllocCPUBuf --> MainLoop
    
    MainLoop[Time Integration Loop] --> RK[Runge-Kutta Scheme:<br/>6 substeps alternating X/Z]
    
    RK --> StepX[Step in X direction]
    StepX --> XTend[xtend: Compute X fluxes]
    
    XTend --> HaloX[exchange_halo_x]
    HaloX --> PackX{OpenACC enabled?}
    
    PackX -->|Yes| PackGPU[!$acc parallel loop:<br/>Pack halo data on GPU]
    PackX -->|No| PackCPU[pack_strip:<br/>Pack halo data on CPU]
    
    PackGPU --> MPISendRecvX
    PackCPU --> MPISendRecvX
    
    MPISendRecvX[MPI_Isend/MPI_Irecv:<br/>Non-blocking send/recv to neighbors]
    MPISendRecvX --> MPIWaitX[MPI_Waitall:<br/>Wait for communication]
    
    MPIWaitX --> UnpackX{OpenACC enabled?}
    UnpackX -->|Yes| UnpackGPU[!$acc parallel loop:<br/>Unpack received data on GPU]
    UnpackX -->|No| UnpackCPU[unpack_strip:<br/>Unpack on CPU]
    
    UnpackGPU --> FluxX
    UnpackCPU --> FluxX
    
    FluxX[Compute X fluxes with:<br/>!$acc parallel loop OR !$omp parallel do]
    FluxX --> UpdateX[Update tendency in X]
    
    UpdateX --> StepZ[Step in Z direction]
    StepZ --> ZTend[ztend: Compute Z fluxes]
    
    ZTend --> HaloZ[exchange_halo_z:<br/>Apply boundary conditions]
    HaloZ --> FluxZ{OpenACC enabled?}
    
    FluxZ -->|Yes| FluxZGPU[!$acc parallel loop:<br/>Compute Z fluxes on GPU]
    FluxZ -->|No| FluxZCPU[!$omp parallel do:<br/>Compute Z fluxes on CPU]
    
    FluxZGPU --> UpdateZ
    FluxZCPU --> UpdateZ
    
    UpdateZ[Update state:<br/>s2 = s0 + dt * tendency]
    UpdateZ --> CheckTime{etime < sim_time?}
    
    CheckTime -->|Yes| CheckOutput{Output time?}
    CheckOutput -->|Yes| WriteOutput[write_record:<br/>Gather and write to NetCDF]
    WriteOutput --> CheckReduce
    CheckOutput -->|No| CheckReduce
    
    CheckReduce{Diagnostic time?}
    CheckReduce -->|Yes| ComputeMass[Compute mass/energy locally]
    ComputeMass --> ReduceOp[MPI_Allreduce:<br/>Sum across all ranks]
    ReduceOp --> RK
    CheckReduce -->|No| RK
    
    CheckTime -->|No| Cleanup
    
    Cleanup[Cleanup and Finalization]
    Cleanup --> AccCleanup{OpenACC enabled?}
    AccCleanup -->|Yes| ExitData[!$acc exit data:<br/>Delete GPU data]
    AccCleanup -->|No| DeallocCPU
    
    ExitData --> DeallocGPU[Deallocate GPU buffers]
    DeallocGPU --> CloseMPI
    DeallocCPU[Deallocate CPU buffers] --> CloseMPI
    
    CloseMPI[MPI_Finalize] --> End([Program End])

    %% Color definitions
    style InitMPI fill:#e1f5ff,color:#000
    style CreateCart fill:#e1f5ff,color:#000
    style MPISendRecvX fill:#e1f5ff,color:#000
    style MPIWaitX fill:#e1f5ff,color:#000
    style ReduceOp fill:#e1f5ff,color:#000
    style CloseMPI fill:#e1f5ff,color:#000

    style DetectGPU fill:#ffe1e1,color:#000
    style AssignGPU fill:#ffe1e1,color:#000
    style CopyToGPU fill:#ffe1e1,color:#000
    style PackGPU fill:#ffe1e1,color:#000
    style UnpackGPU fill:#ffe1e1,color:#000
    style FluxZGPU fill:#ffe1e1,color:#000
    style ExitData fill:#ffe1e1,color:#000
```

