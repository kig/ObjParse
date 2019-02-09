# ObjParse
Fast OBJ parser for WebGL (~150 MB/s on a 871k triangle [dragon.obj](https://casual-effects.com/g3d/data10/index.html#mesh10))

Supports only vertices and normals at this point. Returns flat Float32Arrays that can be used with drawArrays.

## Usage

Load an OBJ file and parse it after the fetch completes.

    const { vertices, normals } = await ObjParse.load(url);
    
There's also a less tested streaming parser that's ~30-50% slower, but should be more interactive-friendly.

Note that there's still a possibly significant blocking computation at the end of the model load to turn the
indexed vertex arrays into flat vertex arrays.

    const { vertices, normals } = await ObjParse.loadStream(url);
    
If you want to get the index arrays, you can do something like this:

    const res = await fetch(url);
    const buf = await res.arrayBuffer();
    const u8 = new Uint8Array(buf);
    
    const {
        vertCount, uvCount, normCount, faceCount, 
        verts, uvs, norms, vidx, uidx, nidx
    } = ObjParse.parse(u8);
    
Note that these arrays have padding at the end, you need to use the <code>*Count</code> variables to iterate through them.
    
    for (let i = 0; i < vertCount; i++) {
        const x = verts[i * 3];
        const y = verts[i * 3 + 1];
        const z = verts[i * 3 + 2];
        ...
    }

## API

    interface ObjParse.Mesh {
        vertices: Float32Array,
        normals: Float32Array,
        uvs: Float32Array
    }
    
    interface ObjParse.IndexedMesh {
        vertCount: number,
        uvCount: number,
        normCount: number,
        faceCount: number,
        verts: Float32Array,
        uvs: Float32Array,
        norms: Float32Array,
        vidx: Uint32Array,
        uidx: Uint32Array,
        nidx: Uint32Array
    }
    
    async ObjParse.load(url:string, computeNormals:boolean = true) : ObjParse.Mesh
    async ObjParse.loadStream(url:string, computeNormals:boolean = true) : ObjParse.Mesh
    
    ObjParse.parse(u8:Uint8Array) : ObjParse.IndexedMesh
    ObjParse.toFlat(mesh:ObjParse.IndexedMesh, computeNormals:boolean = false) : ObjParse.Mesh
    ObjParse.computeNormals(vertices:Float32Array, normals:Float32Array) : void
    
