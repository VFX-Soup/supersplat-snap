import { Container, Label, VectorInput } from '@playcanvas/pcui';
import { Mat4, Quat, Vec3 } from 'playcanvas';

import { EntityTransformOp } from '../edit-ops';
import { Events } from '../events';
import { Scene } from '../scene';
import { Splat } from '../splat';
import { Transform } from '../transform';

// ── Umeyama optimal similarity transform ─────────────────────────────────────
// Finds R, t, s so that dst ≈ s·R·src + t  (least-squares over N point pairs)

type M3 = number[][];

const tr3  = (m: M3): M3 => [[m[0][0],m[1][0],m[2][0]],[m[0][1],m[1][1],m[2][1]],[m[0][2],m[1][2],m[2][2]]];
const mm3  = (A: M3, B: M3): M3 => { const C: M3=[[0,0,0],[0,0,0],[0,0,0]]; for(let i=0;i<3;i++) for(let j=0;j<3;j++) for(let k=0;k<3;k++) C[i][j]+=A[i][k]*B[k][j]; return C; };
const mv3  = (M: M3, v: number[]) => [M[0][0]*v[0]+M[0][1]*v[1]+M[0][2]*v[2], M[1][0]*v[0]+M[1][1]*v[1]+M[1][2]*v[2], M[2][0]*v[0]+M[2][1]*v[1]+M[2][2]*v[2]];
const det3 = (M: M3) => M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);

function svd3(M: M3): { U: M3; S: number[]; V: M3 } {
    const A: M3 = mm3(tr3(M), M).map(r => [...r]) as M3;
    let V: M3 = [[1,0,0],[0,1,0],[0,0,1]];
    for (let iter = 0; iter < 64; iter++) {
        let mx = 0, p = 0, q = 1;
        for (let i = 0; i < 3; i++) for (let j = i+1; j < 3; j++) if (Math.abs(A[i][j]) > mx) { mx = Math.abs(A[i][j]); p=i; q=j; }
        if (mx < 1e-12) break;
        const th = (A[q][q]-A[p][p])/(2*A[p][q]);
        const tv = (th>=0?1:-1)/(Math.abs(th)+Math.sqrt(1+th*th));
        const c = 1/Math.sqrt(1+tv*tv), s2 = tv*c;
        const [App,Aqq,Apq] = [A[p][p],A[q][q],A[p][q]];
        A[p][p]=c*c*App-2*s2*c*Apq+s2*s2*Aqq; A[q][q]=s2*s2*App+2*s2*c*Apq+c*c*Aqq; A[p][q]=A[q][p]=0;
        for (let i=0;i<3;i++) if(i!==p&&i!==q){const[ai,bi]=[A[i][p],A[i][q]];A[i][p]=A[p][i]=c*ai-s2*bi;A[i][q]=A[q][i]=s2*ai+c*bi;}
        for (let i=0;i<3;i++){const[vi,wi]=[V[i][p],V[i][q]];V[i][p]=c*vi-s2*wi;V[i][q]=s2*vi+c*wi;}
    }
    const S = [Math.sqrt(Math.max(0,A[0][0])),Math.sqrt(Math.max(0,A[1][1])),Math.sqrt(Math.max(0,A[2][2]))];
    const MV = mm3(M,V);
    const U: M3 = MV.map((row,_i) => row.map((v,k) => S[k]>1e-12 ? v/S[k] : 0)) as M3;
    return { U, S, V };
}

function umeyama(src: Vec3[], dst: Vec3[]): { R: M3; t: Vec3; s: number; degenerate: boolean } | null {
    const n = Math.min(src.length, dst.length);
    if (n < 3) return null;
    // centroids
    const ms = new Vec3(), md = new Vec3();
    for (let i=0;i<n;i++){ms.add(src[i]);md.add(dst[i]);}
    ms.mulScalar(1/n); md.mulScalar(1/n);
    // centred coords
    const A = src.slice(0,n).map(p=>[p.x-ms.x,p.y-ms.y,p.z-ms.z]);
    const B = dst.slice(0,n).map(p=>[p.x-md.x,p.y-md.y,p.z-md.z]);
    // variance of src
    let vs = 0; for(let i=0;i<n;i++) vs+=A[i][0]**2+A[i][1]**2+A[i][2]**2; vs/=n;
    // cross-covariance H = (1/n) Bᵀ A
    const H: M3=[[0,0,0],[0,0,0],[0,0,0]];
    for(let i=0;i<n;i++) for(let r=0;r<3;r++) for(let c=0;c<3;c++) H[r][c]+=B[i][r]*A[i][c];
    for(let r=0;r<3;r++) for(let c=0;c<3;c++) H[r][c]/=n;
    const {U,S:sv,V} = svd3(H);
    // Proper Umeyama reflection handling: S = diag(1,1, sign(det(U)·det(V)))
    const detU = det3(U), detV = det3(V);
    const ds = [1, 1, (detU * detV < 0) ? -1 : 1];
    const R = mm3(U.map((row,i)=>row.map(v=>v*ds[i])) as M3, tr3(V));
    const s = vs>1e-10 ? (sv[0]*ds[0]+sv[1]*ds[1]+sv[2]*ds[2])/vs : 1;
    const Rm = mv3(R,[ms.x,ms.y,ms.z]);
    // Degenerate if smallest singular value is tiny relative to largest
    // (happens when points are coplanar/collinear → rotation underdetermined)
    const maxSv = Math.max(sv[0], sv[1], sv[2]);
    const minSv = Math.min(sv[0], sv[1], sv[2]);
    const degenerate = maxSv < 1e-9 || (minSv / maxSv) < 1e-4;
    return { R, t: new Vec3(md.x-s*Rm[0], md.y-s*Rm[1], md.z-s*Rm[2]), s, degenerate };
}

// ── Constants ─────────────────────────────────────────────────────────────────

// Distinct colours per pair index (same colour for A1/B1/C1, different for A2/B2/C2 etc)
const PAIR_COLORS = ['#ff4444','#44aaff','#44dd77','#ffcc22','#cc44ff','#ff8844','#44ffee','#ff44aa','#aaff44','#ff9944'];
// Splat label prefix colours for the panel header
const SPLAT_COLORS = ['#ff8888','#88ccff','#88ffaa','#ffdd88','#dd88ff'];

// ── AlignTool ─────────────────────────────────────────────────────────────────

interface AlignPoint {
    splat: Splat;
    worldPos: Vec3;   // live world-space position (editable by user)
}

class AlignTool {
    activate: () => void;
    deactivate: () => void;

    constructor(events: Events, scene: Scene, parent: HTMLElement, canvasContainer: Container) {

        // ── SVG overlay ───────────────────────────────────────────────────
        const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        svg.id = 'align-tool-svg';
        svg.style.cssText = 'position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none;display:none;';
        canvasContainer.dom.appendChild(svg);

        // ── Panel ─────────────────────────────────────────────────────────
        const panel = document.createElement('div');
        panel.id = 'align-tool-panel';
        panel.style.cssText = [
            'position:absolute', 'top:52px', 'right:24px', 'width:270px',
            'background:#1c1c1c', 'border-radius:6px', 'display:none',
            'flex-direction:column', 'z-index:200',
            'box-shadow:0 4px 24px rgba(0,0,0,0.6)',
            'pointer-events:auto', 'overflow:hidden',
            'font-size:12px', 'color:#ccc'
        ].join(';');
        canvasContainer.dom.appendChild(panel);

        // Header
        const hdr = document.createElement('div');
        hdr.style.cssText = 'background:#111;padding:8px 12px;font-weight:bold;font-size:11px;letter-spacing:1px;color:#f60;text-transform:uppercase;';
        hdr.textContent = '⊕ Splat Align';
        panel.appendChild(hdr);

        const body = document.createElement('div');
        body.style.cssText = 'padding:8px 10px;display:flex;flex-direction:column;gap:6px;';
        panel.appendChild(body);

        // Hint
        const hint = document.createElement('div');
        hint.style.cssText = 'font-size:10px;color:#666;line-height:1.5;';
        hint.textContent = 'Select a splat in Scene Manager, then click on it to place points. Match order: A1↔B1, A2↔B2… All splats align to A.';
        body.appendChild(hint);

        // Status
        const status = document.createElement('div');
        status.style.cssText = 'font-size:11px;color:#aaa;background:#111;border-radius:4px;padding:5px 8px;white-space:pre;line-height:1.6;';
        status.textContent = 'No points placed yet.';
        body.appendChild(status);

        // Points list
        const ptListLabel = document.createElement('div');
        ptListLabel.style.cssText = 'font-size:10px;color:#555;text-transform:uppercase;letter-spacing:1px;';
        ptListLabel.textContent = 'Control Points';
        body.appendChild(ptListLabel);

        const ptList = document.createElement('div');
        ptList.style.cssText = 'display:flex;flex-direction:column;gap:2px;max-height:320px;overflow-y:auto;';
        body.appendChild(ptList);

        // Buttons
        const btnRow = document.createElement('div');
        btnRow.style.cssText = 'display:flex;gap:6px;margin-top:2px;';

        const alignBtn = document.createElement('button');
        alignBtn.textContent = 'Align';
        alignBtn.disabled = true;
        alignBtn.style.cssText = 'flex:2;padding:7px;background:#1a5cb0;color:#fff;border:none;border-radius:4px;cursor:pointer;font-size:12px;font-weight:bold;opacity:0.4;';

        const undoBtn = document.createElement('button');
        undoBtn.textContent = 'Undo Last';
        undoBtn.style.cssText = 'flex:1;padding:7px;background:#333;color:#ccc;border:none;border-radius:4px;cursor:pointer;font-size:11px;';

        const clearBtn = document.createElement('button');
        clearBtn.textContent = 'Clear';
        clearBtn.style.cssText = 'padding:7px 10px;background:#333;color:#ccc;border:none;border-radius:4px;cursor:pointer;font-size:11px;';

        btnRow.appendChild(alignBtn);
        btnRow.appendChild(undoBtn);
        btnRow.appendChild(clearBtn);
        body.appendChild(btnRow);

        const resultEl = document.createElement('div');
        resultEl.style.cssText = 'font-size:10px;color:#8d8;line-height:1.6;background:#111;border-radius:4px;padding:5px 8px;white-space:pre;display:none;';
        body.appendChild(resultEl);

        // ── State ─────────────────────────────────────────────────────────
        let isActive = false;
        let points: AlignPoint[] = [];
        let currentSplat: Splat | null = null;

        // ── Helpers ───────────────────────────────────────────────────────

        const splatsInOrder = (): Splat[] => {
            const seen = new Set<Splat>(); const out: Splat[] = [];
            points.forEach(p => { if (!seen.has(p.splat)) { seen.add(p.splat); out.push(p.splat); } });
            return out;
        };

        const ptsForSplat = (s: Splat) => points.filter(p => p.splat === s);
        const splatIdx    = (s: Splat) => splatsInOrder().indexOf(s);
        const pairColor   = (pairIdx: number) => PAIR_COLORS[pairIdx % PAIR_COLORS.length];

        const splatLabel = (s: Splat): string => {
            const si = splatIdx(s);
            return String.fromCharCode(65 + si); // A, B, C …
        };

        // World → normalised screen (0-1)
        const worldToScreen = (wp: Vec3) => {
            const v = wp.clone();
            scene.camera.worldToScreen(v, v);
            return { x: v.x * canvasContainer.dom.clientWidth, y: v.y * canvasContainer.dom.clientHeight };
        };

        // ── Draggable numeric input ───────────────────────────────────────
        // Matches SuperSplat's style: click-drag left/right changes value

        const makeDragInput = (
            axis: string,
            initialVal: number,
            color: string,
            onChange: (v: number) => void
        ): HTMLElement => {
            const wrap = document.createElement('div');
            wrap.style.cssText = 'display:flex;align-items:center;flex:1;gap:1px;min-width:0;';

            const lbl = document.createElement('span');
            lbl.textContent = axis;
            lbl.style.cssText = `font-size:9px;color:${color};width:9px;flex-shrink:0;`;
            wrap.appendChild(lbl);

            const inp = document.createElement('input');
            inp.type = 'number';
            inp.step = '0.0001';
            inp.value = initialVal.toFixed(4);
            inp.style.cssText = [
                'width:100%', 'min-width:0', 'background:#111', 'border:1px solid #333',
                'border-radius:3px', 'color:#ccc', 'font-size:10px', 'padding:2px 3px',
                'cursor:ew-resize', '-moz-appearance:textfield'
            ].join(';');

            // Drag to scrub value
            let dragStart = 0, startVal = 0, dragging = false;
            inp.addEventListener('pointerdown', (e: PointerEvent) => {
                if (e.button !== 0) return;
                dragStart = e.clientX; startVal = parseFloat(inp.value) || 0;
                dragging = false;
                inp.setPointerCapture(e.pointerId);
                e.stopPropagation();
            });
            inp.addEventListener('pointermove', (e: PointerEvent) => {
                if (!(e.buttons & 1)) return;
                const dx = e.clientX - dragStart;
                if (Math.abs(dx) > 2) dragging = true;
                if (!dragging) return;
                const speed = e.shiftKey ? 0.0001 : e.altKey ? 0.01 : 0.001;
                const newVal = startVal + dx * speed;
                inp.value = newVal.toFixed(4);
                onChange(newVal);
                e.stopPropagation();
            });
            inp.addEventListener('pointerup', (e: PointerEvent) => {
                if (dragging) { e.preventDefault(); inp.blur(); }
            });
            inp.addEventListener('change', () => {
                onChange(parseFloat(inp.value) || 0);
            });
            inp.addEventListener('click', (e) => {
                if (dragging) { e.preventDefault(); return; }
                inp.style.cursor = 'text';
            });
            inp.addEventListener('blur', () => { inp.style.cursor = 'ew-resize'; dragging = false; });

            wrap.appendChild(inp);
            wrap.dataset.inp = '1'; // mark for later updates
            (wrap as any)._inp = inp;
            return wrap;
        };

        // ── Points list UI ────────────────────────────────────────────────

        const rebuildPointsList = () => {
            ptList.innerHTML = '';
            const splats = splatsInOrder();
            if (points.length === 0) {
                const empty = document.createElement('div');
                empty.style.cssText = 'font-size:11px;color:#444;padding:4px 2px;';
                empty.textContent = '— no points yet —';
                ptList.appendChild(empty);
                return;
            }

            splats.forEach((splat) => {
                const si = splatIdx(splat);
                const label = splatLabel(splat);
                const scolor = SPLAT_COLORS[si % SPLAT_COLORS.length];

                // Splat header
                const sHdr = document.createElement('div');
                sHdr.style.cssText = `color:${scolor};font-size:10px;text-transform:uppercase;letter-spacing:0.5px;padding:5px 2px 2px;border-top:1px solid #2a2a2a;margin-top:2px;`;
                const splatName = (splat as any)?.name ?? `Splat ${label}`;
                sHdr.textContent = `${label} — ${splatName}${si === 0 ? ' (base)' : ''}`;
                ptList.appendChild(sHdr);

                ptsForSplat(splat).forEach((pt, pi) => {
                    const row = document.createElement('div');
                    row.style.cssText = 'display:flex;align-items:center;gap:3px;padding:2px 0;';

                    // Coloured badge — same colour for same pair index across splats
                    const badge = document.createElement('span');
                    badge.style.cssText = `background:${pairColor(pi)};color:#000;border-radius:50%;width:15px;height:15px;display:flex;align-items:center;justify-content:center;font-size:8px;font-weight:bold;flex-shrink:0;`;
                    badge.textContent = `${label}${pi + 1}`;
                    row.appendChild(badge);

                    // X Y Z drag inputs
                    const axes: Array<[string,string,'x'|'y'|'z']> = [['X','#f88','x'],['Y','#8f8','y'],['Z','#88f','z']];
                    axes.forEach(([axLabel, axColor, axKey]) => {
                        const w = makeDragInput(axLabel, pt.worldPos[axKey], axColor, (val) => {
                            pt.worldPos[axKey] = val;
                            scene.forceRender = true;
                        });
                        row.appendChild(w);
                    });

                    // Delete
                    const del = document.createElement('button');
                    del.textContent = '×';
                    del.style.cssText = 'background:none;border:none;color:#555;font-size:13px;cursor:pointer;padding:0 2px;flex-shrink:0;';
                    del.addEventListener('click', () => {
                        points.splice(points.indexOf(pt), 1);
                        updateAll();
                        scene.forceRender = true;
                    });
                    row.appendChild(del);
                    ptList.appendChild(row);
                });
            });
        };

        const updateStatus = () => {
            const splats = splatsInOrder();
            if (points.length === 0) {
                status.textContent = 'No points placed yet.';
                alignBtn.disabled = true; alignBtn.style.opacity = '0.4';
                return;
            }
            const lines: string[] = [];
            let canAlign = splats.length >= 2;
            splats.forEach((s, i) => {
                const n = ptsForSplat(s).length;
                const name = (s as any)?.name ?? `Splat ${String.fromCharCode(65+i)}`;
                lines.push(`${String.fromCharCode(65+i)} (${name}): ${n} pt${n!==1?'s':''}`);
                if (n < 3) canAlign = false;
            });
            if (splats.length < 2) { lines.push('⚠ Need points on 2 splats'); canAlign = false; }
            else if (!canAlign) { lines.push('⚠ Need ≥3 points per splat'); }
            else { lines.push(`✓ ${Math.min(...splats.map(s=>ptsForSplat(s).length))} pairs — ready`); }
            status.textContent = lines.join('\n');
            alignBtn.disabled = !canAlign;
            alignBtn.style.opacity = canAlign ? '1' : '0.4';
        };

        const updateAll = () => { updateStatus(); rebuildPointsList(); };

        // ── SVG dots & lines ──────────────────────────────────────────────

        const redrawSvg = () => {
            while (svg.lastChild) svg.removeChild(svg.lastChild);
            if (!isActive || points.length === 0) return;

            const splats = splatsInOrder();

            // Dashed lines connecting matched pairs across splats
            if (splats.length >= 2) {
                const baseArr = ptsForSplat(splats[0]);
                for (let si = 1; si < splats.length; si++) {
                    const otherArr = ptsForSplat(splats[si]);
                    const np = Math.min(baseArr.length, otherArr.length);
                    for (let i = 0; i < np; i++) {
                        const a = worldToScreen(baseArr[i].worldPos);
                        const b = worldToScreen(otherArr[i].worldPos);
                        const line = document.createElementNS(svg.namespaceURI, 'line') as SVGLineElement;
                        line.setAttribute('x1',`${a.x}`); line.setAttribute('y1',`${a.y}`);
                        line.setAttribute('x2',`${b.x}`); line.setAttribute('y2',`${b.y}`);
                        line.setAttribute('stroke', pairColor(i));
                        line.setAttribute('stroke-opacity','0.5');
                        line.setAttribute('stroke-width','1.5');
                        line.setAttribute('stroke-dasharray','5,3');
                        svg.appendChild(line);
                    }
                }
            }

            // Dots with labels
            points.forEach(pt => {
                const si = splatIdx(pt.splat);
                const pi = ptsForSplat(pt.splat).indexOf(pt);
                const color = pairColor(pi);
                const s = worldToScreen(pt.worldPos);
                const label = `${String.fromCharCode(65+si)}${pi+1}`;

                const circle = document.createElementNS(svg.namespaceURI,'circle') as SVGCircleElement;
                circle.setAttribute('cx',`${s.x}`); circle.setAttribute('cy',`${s.y}`);
                circle.setAttribute('r','8'); circle.setAttribute('fill',color);
                circle.setAttribute('stroke',SPLAT_COLORS[si%SPLAT_COLORS.length]);
                circle.setAttribute('stroke-width','2');
                svg.appendChild(circle);

                const txt = document.createElementNS(svg.namespaceURI,'text') as SVGTextElement;
                txt.setAttribute('x',`${s.x}`); txt.setAttribute('y',`${s.y+4}`);
                txt.setAttribute('text-anchor','middle'); txt.setAttribute('fill','#000');
                txt.setAttribute('font-size','8'); txt.setAttribute('font-weight','bold');
                txt.textContent = label;
                svg.appendChild(txt);
            });
        };

        // ── Picking ───────────────────────────────────────────────────────

        let clicked = false;
        const pointerdown = (e: PointerEvent) => { if (e.button===0) clicked=true; };
        const pointermove = () => { clicked=false; };
        const pointerup   = async (e: PointerEvent) => {
            if (!clicked || e.button!==0) return;
            clicked = false;
            if (!currentSplat) { status.textContent='Select a splat in Scene Manager first.'; return; }

            const cx = e.offsetX / canvasContainer.dom.clientWidth;
            const cy = e.offsetY / canvasContainer.dom.clientHeight;
            const result = await scene.camera.intersectSplat(currentSplat, cx, cy);
            if (!result) { status.textContent='No surface found — click directly on the selected splat.'; return; }

            // worldPos is the direct hit point in world space — no local transform needed
            points.push({ splat: currentSplat, worldPos: result.position.clone() });
            updateAll();
            scene.forceRender = true;
            e.preventDefault(); e.stopPropagation();
        };

        // ── Alignment ─────────────────────────────────────────────────────

        alignBtn.addEventListener('click', () => {
            const splats = splatsInOrder();
            if (splats.length < 2) return;

            const baseArr = ptsForSplat(splats[0]);
            // dst = base world positions (stays fixed)
            const dst = baseArr.map(p => p.worldPos.clone());

            const results: string[] = [];

            for (let si = 1; si < splats.length; si++) {
                const target = splats[si];
                const targetArr = ptsForSplat(target);
                const n = Math.min(dst.length, targetArr.length);
                // src = current world positions of target's points
                const src = targetArr.slice(0,n).map(p => p.worldPos.clone());
                const dstN = dst.slice(0,n);

                const res = umeyama(src, dstN);
                if (!res) { results.push(`${splatLabel(target)}: need ≥3 pairs`); continue; }
                const { R, t, s, degenerate } = res;

                if (degenerate) {
                    results.push(`${splatLabel(target)}: ⚠ points are coplanar/collinear!\nAdd points at DIFFERENT HEIGHTS (not all on the ground) so the 3D rotation can be solved. Aligned anyway, but rotation may be wrong.`);
                }

                // Save old transform (local == world since splats are children of contentRoot)
                const oldPos   = target.entity.getLocalPosition().clone();
                const oldRot   = target.entity.getLocalRotation().clone();
                const oldScale = target.entity.getLocalScale().clone();

                // Build a PURE rotation quaternion from R (no scale baked in).
                // R is orthonormal from Umeyama, so this is safe.
                const alignRot = new Quat();
                const rotMat = new Mat4();
                const rd = rotMat.data;
                rd[ 0]=R[0][0]; rd[ 1]=R[1][0]; rd[ 2]=R[2][0]; rd[ 3]=0;
                rd[ 4]=R[0][1]; rd[ 5]=R[1][1]; rd[ 6]=R[2][1]; rd[ 7]=0;
                rd[ 8]=R[0][2]; rd[ 9]=R[1][2]; rd[10]=R[2][2]; rd[11]=0;
                rd[12]=0;       rd[13]=0;       rd[14]=0;       rd[15]=1;
                alignRot.setFromMat4(rotMat);
                alignRot.normalize();

                // Compose with existing transform.
                // World transform of a point currently is: world = oldRot·(oldScale·local) + oldPos
                // Alignment applies: aligned = s·alignRot·world + t
                // New TRS:
                //   newScale = s · oldScale
                //   newRot   = alignRot · oldRot
                //   newPos   = s·alignRot·oldPos + t
                const newScale = oldScale.clone().mulScalar(s);

                const newRot = new Quat();
                newRot.mul2(alignRot, oldRot);
                newRot.normalize();

                // newPos = s · (alignRot * oldPos) + t
                const rotatedOldPos = alignRot.transformVector(oldPos.clone());
                const newPos = new Vec3(
                    s * rotatedOldPos.x + t.x,
                    s * rotatedOldPos.y + t.y,
                    s * rotatedOldPos.z + t.z
                );

                target.entity.setLocalPosition(newPos);
                target.entity.setLocalRotation(newRot);
                target.entity.setLocalScale(newScale);

                // Update this splat's stored point worldPos: aligned = s·alignRot·world + t
                ptsForSplat(target).forEach(pt => {
                    const rotated = alignRot.transformVector(pt.worldPos.clone());
                    pt.worldPos.set(
                        s * rotated.x + t.x,
                        s * rotated.y + t.y,
                        s * rotated.z + t.z
                    );
                });

                // Residual error
                let err = 0;
                for (let i=0;i<n;i++) {
                    err += ptsForSplat(target)[i].worldPos.distance(dstN[i]);
                }

                const op = new EntityTransformOp({
                    splat: target,
                    oldt: new Transform(oldPos, oldRot, oldScale),
                    newt: new Transform(newPos, newRot, newScale)
                });
                events.fire('edit.add', op);

                const name = (target as any)?.name ?? `Splat ${splatLabel(target)}`;
                results.push(`${splatLabel(target)} (${name}) aligned ✓\nscale: ${s.toFixed(4)}  mean err: ${(err/n).toFixed(4)}m`);
            }

            resultEl.style.display = 'block';
            resultEl.textContent = results.join('\n') || 'Done.';
            updateAll();
            scene.forceRender = true;
        });

        undoBtn.addEventListener('click', () => {
            if (points.length > 0) { points.pop(); updateAll(); scene.forceRender = true; }
        });

        clearBtn.addEventListener('click', () => {
            points = []; resultEl.style.display='none'; updateAll(); scene.forceRender = true;
        });

        panel.addEventListener('pointerdown', e => e.stopPropagation());

        // ── Events ────────────────────────────────────────────────────────

        events.on('selection.changed', (splat: Splat) => {
            currentSplat = splat?.visible ? splat : null;
        });

        events.on('postrender', () => {
            if (isActive) redrawSvg();
        });

        // ── Activate / Deactivate ─────────────────────────────────────────

        const transformHandler = { activate(){}, deactivate(){} };

        this.activate = () => {
            isActive = true;
            panel.style.display = 'flex';
            svg.style.display = 'block';
            parent.style.display = 'block';
            canvasContainer.dom.addEventListener('pointerdown', pointerdown);
            canvasContainer.dom.addEventListener('pointermove', pointermove);
            canvasContainer.dom.addEventListener('pointerup', pointerup, true);
            updateAll();
            events.fire('transformHandler.push', transformHandler);
        };

        this.deactivate = () => {
            isActive = false;
            panel.style.display = 'none';
            svg.style.display = 'none';
            parent.style.display = 'none';
            canvasContainer.dom.removeEventListener('pointerdown', pointerdown);
            canvasContainer.dom.removeEventListener('pointermove', pointermove);
            canvasContainer.dom.removeEventListener('pointerup', pointerup, true);
            events.fire('transformHandler.pop');
        };
    }
}

export { AlignTool };
