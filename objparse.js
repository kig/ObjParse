/*
   Copyright 2018 Heichen Ltd
   
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

"use strict";

const ObjParse = (function() {
    const MINUS = '-'.charCodeAt(0);
    const DOT = '.'.charCodeAt(0);
    const ZERO = '0'.charCodeAt(0);
    const SPACE = ' '.charCodeAt(0);
    const LF = '\n'.charCodeAt(0);
    const CR = '\r'.charCodeAt(0);
    const SLASH = '/'.charCodeAt(0);
    const V = 'v'.charCodeAt(0);
    const N = 'n'.charCodeAt(0);
    const F = 'f'.charCodeAt(0);
    const E = 'e'.charCodeAt(0);
    const T = 't'.charCodeAt(0);

    const parseVert = function(u8, start, end, verts, vertOff) {
        while (u8[start] === SPACE) {
            start++;
        }
        var num = 0;
        var decimal = 0;
        var decExpMag = 0;
        var sign = 1;
        var exp = 0;
        var expMag = 0;
        var expSign = 1;
        for (let i = start; i < end; ++i) {
            const c = u8[i];
            if (c === SPACE) {
                verts[vertOff] = sign * num * Math.pow(10, decExpMag + expSign * expMag);
                num = 0;
                decimal = 0;
                decExpMag = 0;
                sign = 1;
                exp = 0;
                expMag = 0;
                expSign = 1;
                ++vertOff;
                ++i;
                while (i < end && u8[i] === SPACE) { 
                    ++i;
                }
                --i;
            } else if (c === MINUS) {
                if (exp > 0) {
                    expSign = -1;
                } else {
                    sign = -1;
                }
            } else if (c === DOT) {
                decimal = -1;
            } else if (c === E) {
                exp = 1;
                decimal = 0;
            } else if (exp > 0) {
                expMag = expMag * 10 + (c - ZERO);
            } else {
                decExpMag += decimal;
                num = num * 10 + (c - ZERO);
            }
        }
        verts[vertOff] = sign * num * Math.pow(10, decExpMag + expSign * expMag);
    };

    const parseFace = function(u8, start, end, faceOff, faceTargets, targetCounts) {
        while (u8[start] === 32) {
            start++;
        }
        var num = 0;
        var sign = 1;
        var count = targetCounts[0];
        var target =  faceTargets[0];
        var faceTarget = 0;
        for (let i = start; i < end; ++i) {
            const c = u8[i];
            if (c === SPACE) {
                target[faceOff] = (sign > 0 ? num-1 : count-num);
                num = 0;
                sign = 1;
                faceTarget = 0;
                ++faceOff;
                ++i;
                while (i < end && u8[i] === SPACE) { 
                    ++i;
                }
                --i;
                count = targetCounts[0];
                target = faceTargets[0];
            } else if (c === SLASH) {
                target[faceOff] = (sign > 0 ? num-1 : count-num);
                num = 0;
                sign = 1;
                ++faceTarget;
                ++i;
                while (i < end && u8[i] === SLASH) { 
                    ++i;
                    ++faceTarget;
                }
                --i;
                if (faceTarget > 2) {
                    throw new Error("Only vertex/uv/normal faces supported");
                }
                count = targetCounts[faceTarget];
                target = faceTargets[faceTarget];
            } else if (c === MINUS) {
                sign = -1;
            } else {
                num = num * 10 + (c - ZERO);
            }
        }
        target[faceOff] = (sign > 0 ? num-1 : count-num);
    };

    const parse = function(u8) {
        // Let's round this up to 16k for the heck of it.
        // Most model arrays should fit under this, avoiding reallocs.
        //
        const arrayLen = Math.ceil((u8.byteLength / 20) / (2**12)) * (2**12);
        
        let verts = new Float32Array(arrayLen);
        let uvs = new Float32Array(arrayLen);
        let norms = new Float32Array(arrayLen);
        let vidx = new Uint32Array(arrayLen);
        let uidx = new Uint32Array(arrayLen);
        let nidx = new Uint32Array(arrayLen);

        let vertCount = 0;
        let uvCount = 0;
        let normCount = 0;
        let faceCount = 0;

        const faceTargets = [vidx, uidx, nidx];
        const targetCounts = [vertCount, uvCount, normCount];

        let lineStart = 0;

        for (let i = 0; i < u8.length; ++i) {
            while (i < u8.length && u8[i++] !== LF) { }
                const lineEnd = u8[i-2] === CR ? i-2 : i-1;
                const c0 = u8[lineStart];
                const c1 = u8[lineStart+1];
                const vert = c0 === V;
                const face = c0 === F;
                const norm = vert && c1 === N;
                const uv = vert && c1 === T;
                if (norm) {
                    parseVert(u8, lineStart + 3, lineEnd, norms, normCount*3);
                    ++normCount;
                    if (norms.length <= normCount * 3) {
                        norms = expandFloat32Array(norms);
                    }
                } else if (uv) {
                    parseVert(u8, lineStart + 3, lineEnd, uvs, uvCount*2);
                    ++uvCount;
                    if (uvs.length <= uvCount * 2) {
                        uvs = expandFloat32Array(uvs);
                    }
                } else if (vert) {
                    parseVert(u8, lineStart + 2, lineEnd, verts, vertCount*3);
                    ++vertCount;
                    if (verts.length <= vertCount * 3) {
                        verts = expandFloat32Array(verts);
                    }
                } else if (face) {
                    if  (faceTargets[0] === null) {
                        faceTargets[0] = vidx;
                        faceTargets[1] = uidx;
                        faceTargets[2] = nidx;
                        targetCounts[0] = vertCount;
                        targetCounts[1] = uvCount;
                        targetCounts[2] = normCount;
                    }
                    parseFace(u8, lineStart + 2, lineEnd, faceCount*3, faceTargets, targetCounts);
                    ++faceCount;
                    if (vidx.length <= faceCount * 3) {
                        faceTargets[0] = vidx = expandUint32Array(vidx);
                        faceTargets[1] = nidx = expandUint32Array(nidx);
                        faceTargets[2] = uidx = expandUint32Array(uidx);
                    }
                }
                lineStart = i;
        }
        
        return {
            vertCount, uvCount, normCount, faceCount, 
            verts, uvs, norms, vidx, uidx, nidx
        };
    };

    class ParseState {
        constructor() {
            this.verts = new Float32Array(1048576);
            this.uvs = new Float32Array(1048576);
            this.norms = new Float32Array(1048576);
            this.vidx = new Uint32Array(1048576);
            this.uidx = new Uint32Array(1048576);
            this.nidx = new Uint32Array(1048576);
            this.vertCount = 0;
            this.uvCount = 0;
            this.normCount = 0;
            this.faceCount = 0;
            this.line = new Uint8Array(256);
            this.lineEnd = 0;
            this.lineValid = true;
        }
    }
    
    const parseStream = function(u8, state) {
        for (let offset = 0; offset < u8.length; offset += 4096) {
            const end = offset + Math.min(u8.length - offset, 4096);

            const verts = state.verts;
            const uvs = state.uvs;
            const norms = state.norms;
            const vidx = state.vidx;
            const uidx = state.uidx;
            const nidx = state.nidx;

            let vertCount = state.vertCount;
            let uvCount = state.uvCount;
            let normCount = state.normCount;
            let faceCount = state.faceCount;

            const faceTargets = [vidx, uidx, nidx];
            const targetCounts = [vertCount, uvCount, normCount];

            const line = state.line;
            let lineEnd = state.lineEnd;
            let lineValid = state.lineValid;

            for (let i = offset; i < end; ++i) {
                while (i < end && u8[i] !== LF) {
                    if (lineEnd < line.length) {
                        line[lineEnd++] = u8[i]; 
                    } else {
                        lineValid = false;
                    }
                    ++i;
                }
                if (i < end && u8[i] === LF) {
                    if (lineValid) {
                        const lineEndC = line[lineEnd-2] === CR ? lineEnd-2 : lineEnd-1;
                        const c0 = line[0];
                        const c1 = line[1];
                        const vert = c0 === V;
                        const norm = vert && c1 === N;
                        const uv = vert && c1 === T;
                        const face = c0 === F;
                        if (norm) {
                            parseVert(line, 3, lineEndC, norms, normCount*3);
                            normCount++;
                        } else if (uv) {
                            parseVert(line, 3, lineEndC, uvs, uvCount*2);
                            uvCount++;
                        } else if (vert) {
                            parseVert(line, 2, lineEndC, verts, vertCount*3);
                            vertCount++;
                        } else if (face) {
                            if (targetCounts[0] === 0) {
                                targetCounts[0] = vertCount;
                                targetCounts[1] = uvCount;
                                targetCounts[2] = normCount;
                            }
                            parseFace(line, 2, lineEndC, faceCount*3, faceTargets, targetCounts);
                            faceCount++;
                        }
                    }
                    lineEnd = 0;
                    lineValid = true;
                }
            }

            state.lineEnd = lineEnd;
            state.lineValid = lineValid;
            state.normCount = normCount;
            state.faceCount = faceCount;
            state.vertCount = vertCount;
            state.uvCount = uvCount;

            if (vertCount > verts.length / 3 - 32768) {
                state.verts = expandFloat32Array(verts);
            }
            if (normCount > norms.length / 3 - 32768) {
                state.norms = expandFloat32Array(norms);
            }
            if (uvCount > uvs.length / 2 - 32768) {
                state.uvs = expandFloat32Array(uvs);
            }
            if (faceCount > vidx.length / 3 - 32768) {
                state.vidx = expandUint32Array(vidx);
                state.uidx = expandUint32Array(uidx);
                state.nidx = expandUint32Array(nidx);
            }
        }
    };

    const expandFloat32Array = function(arr) {
        const newArr = new Float32Array(arr.length * 2);
        newArr.set(arr);
        return newArr;
    };

    const expandUint32Array = function(arr) {
        const newArr = new Uint32Array(arr.length * 2);
        newArr.set(arr);
        return newArr;
    };

    const computeNormals = function(v, n) {
        for (let j = 0; j < v.length; j += 9) {
            let i = j;
            const x0 = v[i];
            const y0 = v[++i];
            const z0 = v[++i]; 
            const x1 = v[++i];
            const y1 = v[++i];
            const z1 = v[++i]; 
            const x2 = v[++i];
            const y2 = v[++i];
            const z2 = v[++i];
            const e0x = x1-x0;
            const e0y = y1-y0;
            const e0z = z1-z0;
            const e1x = x2-x0;
            const e1y = y2-y0;
            const e1z = z2-z0;
            const nx = e0y*e1z - e0z*e1y;
            const ny = e0z*e1x - e0x*e1z;
            const nz = e0x*e1y - e0y*e1x;
            const idn = 1.0 / Math.sqrt(nx*nx + ny*ny + nz*nz);
            const nnx = nx * idn;
            const nny = ny * idn;
            const nnz = nz * idn;
            i = j;
            n[i] = nnx;
            n[++i] = nny;
            n[++i] = nnz;
            n[++i] = nnx;
            n[++i] = nny;
            n[++i] = nnz;
            n[++i] = nnx;
            n[++i] = nny;
            n[++i] = nnz;
        }
    };

    const _toFlatVerts = function(obj, normals) {
        const verts = obj.verts;
        const vidx = obj.vidx;
        const faceCount = obj.faceCount;
        const faceCount3 = faceCount * 3;
        const flatVerts = new Float32Array( 3 * faceCount3 );
        const flatNorms = new Float32Array( 3 * faceCount3 );
        const flatUVs = new Float32Array( 2 * faceCount3 );
        for (let i = 0, j = 0; i < faceCount3; ++i, j += 3) {
            const vi = vidx[i] * 3;
            flatVerts[j] = verts[vi];
            flatVerts[j+1] = verts[vi+1];
            flatVerts[j+2] = verts[vi+2];
        }
        if (normals) {
            console.time('ObjParse.computeNormals');
            computeNormals(flatVerts, flatNorms);
            console.timeEnd('ObjParse.computeNormals');
        }
        return { vertices: flatVerts, uvs: flatUVs, normals: flatNorms };
    };

    const _toFlatVertsNorms = function(obj) {
        const verts = obj.verts;
        const norms = obj.norms;
        const vidx = obj.vidx;
        const nidx = obj.nidx;
        const faceCount = obj.faceCount;
        const faceCount3 = faceCount * 3;
        const flatVerts = new Float32Array( 3 * faceCount3 );
        const flatNorms = new Float32Array( 3 * faceCount3 );
        const flatUVs = new Float32Array( 2 * faceCount3 );
        for (let i = 0, j = 0; i < faceCount3; ++i, j += 3) {
            const vi = vidx[i] * 3;
            const ni = nidx[i] * 3;
            flatVerts[j] = verts[vi];
            flatVerts[j+1] = verts[vi+1];
            flatVerts[j+2] = verts[vi+2];
            flatNorms[j] = norms[ni];
            flatNorms[j+1] = norms[ni+1];
            flatNorms[j+2] = norms[ni+2];
        }
        return { vertices: flatVerts, uvs: flatUVs, normals: flatNorms };
    };

    const _toFlatVertsUVs = function(obj, normals) {
        const verts = obj.verts;
        const uvs = obj.uvs;
        const vidx = obj.vidx;
        const uidx = obj.uidx;
        const faceCount = obj.faceCount;
        const faceCount3 = faceCount * 3;
        const flatVerts = new Float32Array( 3 * faceCount3 );
        const flatNorms = new Float32Array( 3 * faceCount3 );
        const flatUVs = new Float32Array( 2 * faceCount3 );
        for (let i = 0, j = 0, k = 0; i < faceCount3; ++i, j += 3, k += 2) {
            const vi = vidx[i] * 3;
            const ui = uidx[i] * 2;
            flatVerts[j] = verts[vi];
            flatVerts[j+1] = verts[vi+1];
            flatVerts[j+2] = verts[vi+2];
            flatUVs[k] = uvs[ui];
            flatUVs[k+1] = uvs[ui+1];
        }
        if (normals) {
            console.time('ObjParse.computeNormals');
            computeNormals(flatVerts, flatNorms);
            console.timeEnd('ObjParse.computeNormals');
        }
        return { vertices: flatVerts, uvs: flatUVs, normals: flatNorms };
    };

    const _toFlatVertsNormsUVs = function(obj) {
        const verts = obj.verts;
        const uvs = obj.uvs;
        const norms = obj.norms;
        const vidx = obj.vidx;
        const uidx = obj.uidx;
        const nidx = obj.nidx;
        const faceCount = obj.faceCount;
        const faceCount3 = faceCount * 3;
        const flatVerts = new Float32Array( 3 * faceCount3 );
        const flatNorms = new Float32Array( 3 * faceCount3 );
        const flatUVs = new Float32Array( 2 * faceCount3 );
        for (let i = 0, j = 0, k = 0; i < faceCount3; ++i, j += 3, k += 2) {
            const vi = vidx[i] * 3;
            const ui = uidx[i] * 2;
            const ni = nidx[i] * 3;
            flatVerts[j] = verts[vi];
            flatVerts[j+1] = verts[vi+1];
            flatVerts[j+2] = verts[vi+2];
            flatNorms[j] = norms[ni];
            flatNorms[j+1] = norms[ni+1];
            flatNorms[j+2] = norms[ni+2];
            flatUVs[k] = uvs[ui];
            flatUVs[k+1] = uvs[ui+1];
        }
        return { vertices: flatVerts, uvs: flatUVs, normals: flatNorms };
    };

    const toFlat = function(obj, normals = false) {
        if (obj.normCount === 0 && obj.uvCount === 0) {
            return _toFlatVerts(obj, normals);
        } else if (obj.uvCount === 0) {
            return _toFlatVertsNorms(obj);
        } else if (obj.normCount === 0) {
            return _toFlatVertsUVs(obj, normals);
        } else {
            return _toFlatVertsNormsUVs(obj);
        }
    };

    const loadStream = async function(url, computeNormals = true) {
        console.time('ObjParse.loadStream');

        console.time('ObjParse.fetch');

        const res = await fetch(url);
        const body = await res.body;
        const reader = body.getReader();
        const parseState = new ParseState();

        while (true) {
            const read = await reader.read();
            const done = read.done;
            const value = read.value;

            if (done) {
                console.timeEnd('ObjParse.fetch');
                console.time('ObjParse.toFlat');
                const flat = toFlat(parseState, computeNormals);
                console.timeEnd('ObjParse.toFlat');
                console.timeEnd('ObjParse.loadStream');
                return flat;
            }

            parseStream(value, parseState);

        }
    };

    const load = async function(url, computeNormals = true) {
        console.time('ObjParse.load');

        console.time('ObjParse.fetch');
        const res = await fetch(url);
        const buf = await res.arrayBuffer();
        const u8 = new Uint8Array(buf);
        console.timeEnd('ObjParse.fetch');

        console.time('ObjParse.parse');
        const parseState = parse(u8);
        console.timeEnd('ObjParse.parse');

        console.time('ObjParse.toFlat');
        const flat = toFlat(parseState, computeNormals);
        console.timeEnd('ObjParse.toFlat');

        console.timeEnd('ObjParse.load');
        return flat;
    };

    return { load, loadStream, parse, parseStream, toFlat, computeNormals, ParseState };
})();

module.exports = ObjParse;
