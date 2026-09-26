const API = {
  async request(method, path, body, isFormData = false) {
    const options = { method };
    if (body) {
      if (isFormData) {
        options.body = body;
      } else {
        options.headers = { 'Content-Type': 'application/json' };
        options.body = JSON.stringify(body);
      }
    }
    const response = await fetch(`/api${path}`, options);
    if (!response.ok) {
      const error = await response.json().catch(() => ({ detail: response.statusText }));
      throw new Error(error.detail || 'Request failed');
    }
    return response.json();
  },
  listJobs() {
    return this.request('GET', '/jobs');
  },
  getKnowledgeStats() {
    return this.request('GET', '/knowledge/stats');
  },
  searchKnowledge(q, type, limit = 20) {
    const p = new URLSearchParams({ q: q || '', limit: String(limit) });
    if (type) p.set('type', type);
    return this.request('GET', `/knowledge/search?${p.toString()}`);
  },
  listTrajectories() {
    return this.request('GET', '/trajectories');
  },
  listDesignJobs() {
    return this.request('GET', '/design/jobs');
  },
  createDesignJob(payload) {
    return this.request('POST', '/design/jobs', payload);
  },
  getDesignJob(id) {
    return this.request('GET', `/design/jobs/${id}`);
  },
  getDesignContext(id) {
    return this.request('GET', `/design/jobs/${id}/context`);
  },
  getDesignDrafts(id) {
    return this.request('GET', `/design/jobs/${id}/drafts`);
  },
  deleteDesignDraft(id, draftIdx) {
    return this.request('DELETE', `/design/jobs/${id}/drafts/${draftIdx}`);
  },
  deleteDesignJob(id) {
    return this.request('DELETE', `/design/jobs/${id}`);
  },
  runDesignPipeline(id) {
    return this.request('POST', `/design/jobs/${id}/run_pipeline`);
  },
  getDesignPipelineStatus(id) {
    return this.request('GET', `/design/jobs/${id}/run_pipeline/status`);
  },
  getDesignPipeline(id) {
    return this.request('GET', `/design/jobs/${id}/pipeline`);
  },
  getJob(id) {
    return this.request('GET', `/jobs/${id}`);
  },
  createJob(formData) {
    return this.request('POST', '/jobs', formData, true);
  },
  createCompareJob(payload) {
    return this.request('POST', '/compare/jobs', payload);
  },
  getCompareSummary(id) {
    return this.request('GET', `/compare/jobs/${id}/summary`);
  },
  getCompareTable(id) {
    return this.request('GET', `/compare/jobs/${id}/table`);
  },
  getCompareSubTable(id, name) {
    return this.request('GET', `/compare/jobs/${id}/tables/${name}`);
  },
  listComparePlots(id) {
    return this.request('GET', `/compare/jobs/${id}/plots`);
  },
  startJob(id) {
    return this.request('POST', `/jobs/${id}/start`);
  },
  deleteJob(id) {
    return this.request('DELETE', `/jobs/${id}`);
  },
  getLlmSettings() {
    return this.request('GET', '/settings/llm');
  },
  updateLlmSettings(payload) {
    return this.request('PUT', '/settings/llm', payload);
  },
  testLlmConnection() {
    return this.request('POST', '/settings/llm/test');
  },
  askJob(id, payload) {
    return this.request('POST', `/jobs/${id}/assistant/query`, payload);
  },
  listJobSkills(id) {
    return this.request('GET', `/jobs/${id}/assistant/skills`);
  },
  invokeJobSkill(id, payload) {
    return this.request('POST', `/jobs/${id}/assistant/invoke`, payload);
  },
};

const MODULE_PRESETS = {
  standard: 'quality,identity,contact,rmsf,bsa,rrcs,cluster,landscape,angles,dihedrals,report',
  contact: 'identity,contact,rrcs,report',
  fast: 'quality,identity,contact,report',
  geometry: 'quality,identity,angles,dihedrals,report',
  custom: 'quality,identity,contact,report',
};

const ANALYSIS_MODULES = [
  {
    id: 'quality',
    label: 'Quality',
    description: 'Trajectory integrity, frame counts, and basic validation.',
  },
  {
    id: 'identity',
    label: 'Identity',
    description: 'Chain/component annotation for TCR, peptide, and HLA.',
  },
  {
    id: 'contact',
    label: 'Contacts',
    description: 'Interface contact families and occupancy summaries.',
  },
  {
    id: 'rmsf',
    label: 'RMSF',
    description: 'Residue flexibility and region-level mobility.',
  },
  {
    id: 'bsa',
    label: 'BSA',
    description: 'Interface buried surface area over the trajectory.',
  },
  {
    id: 'rrcs',
    label: 'RRCS',
    description: 'Residue-residue contact score and hotspot ranking.',
  },
  {
    id: 'cluster',
    label: 'Cluster',
    description: 'Interface conformational clustering and dominant states.',
  },
  {
    id: 'landscape',
    label: 'Landscape',
    description: 'Low-dimensional conformational landscape (PCA/UMAP/TICA).',
  },
  {
    id: 'angles',
    label: 'Docking angles',
    description: 'TCR-pMHC binding geometry — crossing / incident / tilt angles.',
  },
  {
    id: 'dihedrals',
    label: 'Dihedrals (φ/ψ)',
    description: 'Per-residue Ramachandran analysis and backbone flexibility ranking.',
  },
  {
    id: 'report',
    label: 'Report',
    description: 'Integrated HTML report and web-readable result index.',
  },
];

const QUICK_QUESTIONS = [
  'Summarize this trajectory',
  'Diagnose stability',
  'Summarize interface',
  'Identify hotspots',
];

const App = {
  agentClient: null,
  agentHistory: null,
  async start() {
    window.addEventListener('hashchange', () => this.navigate());
    await this.navigate();
    window.setInterval(() => {
      const hash = location.hash.slice(1) || '/';
      if (hash === '/jobs' || hash.startsWith('/jobs/')) this.navigate(false);
    }, 6000);
  },
  async navigate(showLoading = true) {
    const fullHash = (location.hash.slice(1) || '/').replace(/\/$/, '') || '/';
    const [hashPath, hashQuery] = fullHash.split('?');
    const hash = hashPath || '/';
    const queryParams = new URLSearchParams(hashQuery || '');
    this.queryParams = queryParams;
    const content = document.getElementById('page-content');
    if (showLoading) content.innerHTML = '<div class="panel">Loading...</div>';
    document.querySelectorAll('.nav a').forEach(link => {
      link.classList.toggle('active', link.getAttribute('href') === `#${hash}`);
    });

    // Save agent history when leaving agent page
    if (this.agentClient && hash !== '/agent') {
      this.agentHistory = this.agentClient.saveHistory();
      this.agentClient.cleanup();
      this.agentClient = null;
    }

    // Clean up ImmunoScope (design) WebSocket when leaving conversation
    if (window._designCleanup && !hash.startsWith('/copilot/')) {
      try { window._designCleanup(); } catch (e) {}
      window._designCleanup = null;
    }

    if (hash === '/agent') {
      document.getElementById('page-title').textContent = 'Query-Agent';
      content.innerHTML = renderAgent();
      this.agentClient = new AgentClient();
      // Restore history if exists
      if (this.agentHistory) {
        this.agentClient.restoreHistory(this.agentHistory);
      }
      return;
    }
    if (hash === '/new') {
      document.getElementById('page-title').textContent = 'New Analysis';
      content.innerHTML = await renderNewAnalysis();
      mountNewAnalysis();
      return;
    }
    if (hash === '/jobs') {
      document.getElementById('page-title').textContent = 'Jobs';
      content.innerHTML = await renderJobs();
      return;
    }
    if (hash === '/projects') {
      document.getElementById('page-title').textContent = 'Projects';
      content.innerHTML = renderComingSoon('Projects', 'Organize TCR-pMHC systems, replicate sets, and design campaigns.');
      return;
    }
    if (hash === '/compare') {
      document.getElementById('page-title').textContent = 'Compare';
      content.innerHTML = await renderCompare();
      mountCompare();
      return;
    }
    if (hash === '/reports') {
      document.getElementById('page-title').textContent = 'Reports';
      content.innerHTML = await renderReports();
      return;
    }
    if (hash === '/trajectories') {
      document.getElementById('page-title').textContent = 'Trajectories';
      content.innerHTML = await renderTrajectories();
      return;
    }
    if (hash === '/copilot') {
      document.getElementById('page-title').textContent = 'ImmunoScope';
      content.innerHTML = await renderDesignCopilot();
      mountDesignCopilot();
      return;
    }
    if (hash.startsWith('/copilot/new')) {
      document.getElementById('page-title').textContent = 'New ImmunoScope Session';
      content.innerHTML = await renderNewDesign();
      mountNewDesign();
      return;
    }
    if (hash.startsWith('/copilot/') && hash.includes('/showcase')) {
      const designId = hash.split('/')[2];
      document.getElementById('page-title').textContent = 'ImmunoScope · Recommendations';
      document.body.classList.add('showcase-mode');
      content.innerHTML = await renderDesignShowcase(designId);
      mountDesignShowcase(designId);
      return;
    }
    // Remove showcase mode when leaving the showcase page
    document.body.classList.remove('showcase-mode');
    if (hash.startsWith('/copilot/')) {
      const designId = hash.split('/')[2];
      document.getElementById('page-title').textContent = 'ImmunoScope';
      content.innerHTML = await renderDesignDetail(designId);
      mountDesignDetail(designId);
      return;
    }
    if (hash === '/knowledge') {
      document.getElementById('page-title').textContent = 'Knowledge Base';
      content.innerHTML = renderKnowledgeBase();
      mountKnowledgeBase();
      return;
    }
    if (hash === '/settings') {
      document.getElementById('page-title').textContent = 'Settings';
      content.innerHTML = await renderSettings();
      mountSettings();
      return;
    }
    if (hash === '/docs') {
      document.getElementById('page-title').textContent = 'Documentation';
      content.innerHTML = renderDocs();
      return;
    }
    if (hash.startsWith('/jobs/')) {
      const id = hash.split('/')[2];
      // Detect compare jobs (profile=sampling_compare) and route to compare view
      const job = await API.getJob(id).catch(() => null);
      const isCompareJob = job && job.config && (
        (job.config.modules || []).includes('compare') ||
        String(job.config.profile || '').toLowerCase().includes('compare')
      );
      if (isCompareJob) {
        document.getElementById('page-title').textContent = 'Comparison';
        content.innerHTML = await renderCompareDetail(id, job);
        mountCompareDetail(id);
        return;
      }
      document.getElementById('page-title').textContent = 'Job Detail';
      content.innerHTML = await renderJobDetail(id);
      mountJobDetail(id);
      return;
    }
    document.getElementById('page-title').textContent = 'Dashboard';
    content.innerHTML = await renderDashboard();
  },
};

function escapeHtml(value) {
  const div = document.createElement('div');
  div.textContent = value == null ? '' : String(value);
  return div.innerHTML;
}

function renderComingSoon(title, description) {
  return `
    <div class="panel">
      <span class="eyebrow">Workspace module</span>
      <h2>${escapeHtml(title)}</h2>
      <p class="muted">${escapeHtml(description)}</p>
    </div>
  `;
}

// ---------- Knowledge Base ----------
// Browsable / searchable view of the extracted TCR literature fact blocks
// (ChromaDB `tcr_facts`, served by /api/knowledge). Semantic search returns
// fact cards with the claim + a verbatim evidence quote + PMID.
const KB_TYPES = {
  mutation_effect:        { label: 'Mutation → effect', color: '#15803d', bg: '#dcfce7' },
  binding_measurement:    { label: 'Binding measurement', color: '#0e7490', bg: '#cffafe' },
  structural_observation: { label: 'Structural', color: '#6d28d9', bg: '#ede9fe' },
  design_heuristic:       { label: 'Design heuristic', color: '#a16207', bg: '#fef3c7' },
};
const KB_EXAMPLES = [
  'aromatic substitution in CDR3 increases affinity',
  'CD8 coreceptor affinity and peptide specificity',
  'KD measured by SPR for TCR-pMHC',
  'thermostability effect of CDR loop mutation',
];

function renderKnowledgeBase() {
  const chips = ['', ...Object.keys(KB_TYPES)].map(t => {
    const label = t === '' ? 'All' : KB_TYPES[t].label;
    return `<button class="kb-chip" data-kb-type="${t}" style="padding:5px 12px;border-radius:999px;border:1px solid var(--line);background:#fff;color:var(--text);font-size:12px;font-weight:700;cursor:pointer;">${escapeHtml(label)} <span class="kb-chip-count" style="color:var(--muted);font-weight:600;"></span></button>`;
  }).join('');
  const examples = KB_EXAMPLES.map(e =>
    `<button class="kb-example" data-kb-q="${escapeHtml(e)}" style="background:none;border:0;color:var(--accent-strong);font-size:12px;cursor:pointer;text-decoration:underline;text-underline-offset:2px;padding:2px 0;text-align:left;">${escapeHtml(e)}</button>`
  ).join('<br>');
  return `
    <div class="panel" style="margin-bottom:14px;">
      <span class="eyebrow">Knowledge Base</span>
      <h2 style="margin:2px 0 4px;">Literature fact library</h2>
      <p class="muted" id="kb-subtitle" style="margin:0;">Citable biochemical facts mined from full-text TCR-pMHC papers — search the claim, see the verbatim quote + PMID.</p>
      <div style="display:flex;gap:8px;margin-top:14px;">
        <input id="kb-q" type="text" placeholder="Search facts — e.g. CDR3 aromatic substitution affinity"
               style="flex:1;padding:10px 14px;border-radius:10px;border:1px solid var(--line);font-size:14px;">
        <button id="kb-go" style="padding:10px 20px;border-radius:10px;border:0;background:#1c473c;color:#fffdf8;font-weight:800;font-size:13px;cursor:pointer;">Search</button>
      </div>
      <div id="kb-chips" style="display:flex;flex-wrap:wrap;gap:6px;margin-top:12px;">${chips}</div>
    </div>
    <div id="kb-results">
      <div class="panel" style="color:var(--muted);">
        <div style="font-size:13px;margin-bottom:8px;">Try an example:</div>
        ${examples}
      </div>
    </div>
  `;
}

function renderKbCard(c) {
  const meta = KB_TYPES[c.block_type] || { label: c.block_type || 'fact', color: '#64748b', bg: '#f1f5f9' };
  const confColor = c.confidence === 'high' ? '#15803d' : (c.confidence === 'low' ? '#b91c1c' : '#a16207');
  const pmid = c.pmid
    ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${escapeHtml(c.pmid)}/" target="_blank" rel="noopener" style="padding:2px 8px;border-radius:6px;background:#a36b00;color:#fffdf8;font-size:10px;font-weight:800;text-decoration:none;">PMID ${escapeHtml(c.pmid)}</a>`
    : '';
  const sys = c.system ? `<span class="muted" style="font-size:11px;">${escapeHtml(c.system)}</span>` : '';
  const reg = c.region ? `<span class="muted" style="font-size:11px;">${escapeHtml(c.region)}</span>` : '';
  return `
    <div style="background:#fff;border:1px solid var(--line);border-radius:14px;padding:14px 16px;margin-bottom:10px;">
      <div style="display:flex;align-items:center;gap:8px;margin-bottom:8px;flex-wrap:wrap;">
        <span style="padding:3px 9px;border-radius:999px;background:${meta.bg};color:${meta.color};font-size:10px;font-weight:800;letter-spacing:0.04em;text-transform:uppercase;">${escapeHtml(meta.label)}</span>
        <span style="font-size:11px;font-weight:700;color:${confColor};text-transform:uppercase;">${escapeHtml(c.confidence || '')}</span>
        ${reg}${sys}
        ${c.score != null ? `<span class="muted" style="margin-left:auto;font-size:10px;font-family:ui-monospace,monospace;">match ${c.score}</span>` : ''}
      </div>
      <div style="font-size:14px;font-weight:600;line-height:1.5;color:var(--text);">${escapeHtml(c.statement || '')}</div>
      ${c.evidence_quote ? `<div style="margin-top:8px;padding:8px 12px;border-left:3px solid ${meta.color};background:#faf8f3;border-radius:0 8px 8px 0;font-size:12.5px;line-height:1.55;color:#44504a;font-style:italic;">“${escapeHtml(c.evidence_quote)}”</div>` : ''}
      <div style="display:flex;align-items:center;gap:8px;margin-top:8px;">
        ${pmid}
        ${c.source_section ? `<span class="muted" style="font-size:11px;">§ ${escapeHtml(c.source_section)}</span>` : ''}
      </div>
    </div>
  `;
}

function mountKnowledgeBase() {
  const input = document.getElementById('kb-q');
  const go = document.getElementById('kb-go');
  const chipsWrap = document.getElementById('kb-chips');
  const results = document.getElementById('kb-results');
  if (!input || !results) return;
  let activeType = '';

  // Fill chip counts + total from stats.
  API.getKnowledgeStats().then(s => {
    const sub = document.getElementById('kb-subtitle');
    if (sub) sub.textContent = `${s.total} citable facts from full-text TCR-pMHC papers — search the claim, see the verbatim quote + PMID.`;
    chipsWrap.querySelectorAll('.kb-chip').forEach(btn => {
      const t = btn.getAttribute('data-kb-type');
      const cnt = t === '' ? s.total : (s.by_type[t] || 0);
      const el = btn.querySelector('.kb-chip-count');
      if (el) el.textContent = `(${cnt})`;
    });
  }).catch(() => {});

  function paintChips() {
    chipsWrap.querySelectorAll('.kb-chip').forEach(btn => {
      const on = btn.getAttribute('data-kb-type') === activeType;
      btn.style.background = on ? '#1c473c' : '#fff';
      btn.style.color = on ? '#fffdf8' : 'var(--text)';
      btn.style.borderColor = on ? '#1c473c' : 'var(--line)';
    });
  }
  paintChips();

  async function doSearch() {
    const q = (input.value || '').trim();
    if (!q) return;
    results.innerHTML = `<div class="panel muted">Searching…</div>`;
    try {
      const data = await API.searchKnowledge(q, activeType, 25);
      const cards = data.results || [];
      if (!cards.length) {
        results.innerHTML = `<div class="panel muted">No facts matched “${escapeHtml(q)}”${activeType ? ' in ' + escapeHtml(KB_TYPES[activeType].label) : ''}.</div>`;
        return;
      }
      results.innerHTML =
        `<div class="muted" style="font-size:12px;margin:0 2px 10px;">${cards.length} facts${activeType ? ' · ' + escapeHtml(KB_TYPES[activeType].label) : ''} · ranked by relevance</div>`
        + cards.map(renderKbCard).join('');
    } catch (e) {
      results.innerHTML = `<div class="panel" style="color:#b91c1c;">Search failed: ${escapeHtml(e.message || String(e))}</div>`;
    }
  }

  go.addEventListener('click', doSearch);
  input.addEventListener('keydown', e => { if (e.key === 'Enter') doSearch(); });
  chipsWrap.addEventListener('click', e => {
    const btn = e.target.closest('.kb-chip');
    if (!btn) return;
    activeType = btn.getAttribute('data-kb-type') || '';
    paintChips();
    if ((input.value || '').trim()) doSearch();
  });
  results.addEventListener('click', e => {
    const ex = e.target.closest('.kb-example');
    if (!ex) return;
    input.value = ex.getAttribute('data-kb-q') || '';
    doSearch();
  });
  input.focus();
}

function statusBadge(status) {
  return `<span class="status ${escapeHtml(status)}">${escapeHtml(status)}</span>`;
}

function formatTime(value) {
  return value ? new Date(value).toLocaleString() : '-';
}

async function renderDashboard() {
  const jobs = await API.listJobs().catch(() => []);
  const completed = jobs.filter(job => job.status === 'completed').length;
  const running = jobs.filter(job => job.status === 'running').length;
  const failed = jobs.filter(job => job.status === 'failed').length;
  return `
    <div class="grid cards">
      <div class="card"><span>Total jobs</span><strong>${jobs.length}</strong></div>
      <div class="card"><span>Running</span><strong>${running}</strong></div>
      <div class="card"><span>Completed</span><strong>${completed}</strong></div>
      <div class="card"><span>Failed</span><strong>${failed}</strong></div>
    </div>
    <div class="panel hero-panel">
      <h2>Product Entry</h2>
      <p class="muted">Create an ImmunoScope job from raw MD inputs, run the CLI workflow, inspect the generated report, and ask the Assistant about evidence-backed findings.</p>
      <div class="actions">
        <a class="button" href="#/new">Create analysis job</a>
        <a class="button secondary" href="#/compare">Compare conditions</a>
        <a class="button secondary" href="#/jobs">View jobs</a>
      </div>
    </div>
  `;
}

async function renderCompare() {
  const jobs = await API.listJobs().catch(() => []);
  // Only completed single-system analysis jobs are eligible as compare inputs
  const eligible = jobs.filter(j => j.status === 'completed' && j.run_dir);

  if (eligible.length < 2) {
    return `
      <div class="panel">
        <h2>Compare Systems</h2>
        <p class="muted">You need at least two completed analysis jobs to run a comparison. Currently: <strong>${eligible.length}</strong> available.</p>
        <a class="button" href="#/new">Create new analysis</a>
      </div>
    `;
  }

  // Pre-list job options
  const jobOption = (j, selected = false) =>
    `<option value="${escapeHtml(j.id)}" ${selected ? 'selected' : ''}>${escapeHtml(j.name || j.id)} <span style="color: var(--muted);">(${escapeHtml(j.id)})</span></option>`;

  const optionsA = eligible.map(j => jobOption(j)).join('');
  const optionsB = eligible.map((j, i) => jobOption(j, i === 1)).join('');

  return `
    <form id="compare-form" class="panel analysis-form">
      <div class="form-header">
        <div>
          <h2>Compare Two Systems</h2>
          <p class="muted">Compare two completed analyses. Use this for Standard MD vs Enhanced sampling, replicate checks, WT vs mutant, or different binder/peptide variants.</p>
        </div>
        <label class="switch-row">
          <input name="auto_start" type="checkbox" checked value="true">
          Start immediately
        </label>
      </div>

      <!-- Source selection: side-by-side dropdowns -->
      <div class="grid two-col" style="margin: 8px 0 16px;">
        <div style="padding: 16px; border-radius: 16px; border: 1px solid var(--line); background: linear-gradient(180deg, rgba(245, 241, 232, 0.55) 0%, transparent 100%);">
          <span style="font-size: 10px; font-weight: 800; color: var(--muted); letter-spacing: 0.12em; text-transform: uppercase;">Case A — reference</span>
          <label style="margin-top: 10px;">
            Source analysis
            <select name="case_a_id" id="case-a-select" required>${optionsA}</select>
          </label>
          <label>
            Label A
            <input name="label_a" value="Condition A" placeholder="e.g. Standard MD">
          </label>
        </div>
        <div style="padding: 16px; border-radius: 16px; border: 1px solid var(--line); background: linear-gradient(180deg, rgba(219, 233, 225, 0.40) 0%, transparent 100%);">
          <span style="font-size: 10px; font-weight: 800; color: var(--accent-strong); letter-spacing: 0.12em; text-transform: uppercase;">Case B — comparison</span>
          <label style="margin-top: 10px;">
            Source analysis
            <select name="case_b_id" id="case-b-select" required>${optionsB}</select>
          </label>
          <label>
            Label B
            <input name="label_b" value="Condition B" placeholder="e.g. Enhanced sampling">
          </label>
        </div>
      </div>

      <label>
        Compare job name (optional)
        <input name="job_name" placeholder="e.g. wt_vs_mutant_comparison">
      </label>

      <div class="grid three-col">
        <label>
          Mode
          <select name="comparison_mode">
            <option value="generic">Generic</option>
            <option value="sampling">Sampling</option>
            <option value="replicate">Replicate</option>
            <option value="mutation">Mutation</option>
          </select>
        </label>
        <label>
          Scope
          <select name="comparison_scope">
            <option value="same-system">Same system</option>
            <option value="cross-system">Cross system</option>
            <option value="auto">Auto</option>
          </select>
        </label>
        <label>
          Alignment selection
          <input name="alignment_selection" value="phla_core_ca">
        </label>
      </div>

      <div class="grid two-col">
        <label>
          Residue mapping policy
          <select name="residue_mapping">
            <option value="auto">Auto</option>
            <option value="identity">Identity (strict)</option>
            <option value="region-only">Region only</option>
          </select>
        </label>
      </div>

      <label>
        Notes / context (optional)
        <textarea name="comparison_context" placeholder="What are you comparing and why?"></textarea>
      </label>

      <button type="submit">Create comparison</button>
      <p id="compare-status" class="muted"></p>
    </form>
  `;
}

function mountCompare() {
  const form = document.getElementById('compare-form');
  if (!form) return;
  const status = document.getElementById('compare-status');

  form.addEventListener('submit', async event => {
    event.preventDefault();
    status.textContent = 'Creating compare job...';

    const data = new FormData(form);
    const caseAId = data.get('case_a_id');
    const caseBId = data.get('case_b_id');

    if (caseAId === caseBId) {
      status.textContent = 'Pick two different source analyses for Case A and Case B.';
      return;
    }

    // Resolve job IDs → run_dir paths
    try {
      const [jobA, jobB] = await Promise.all([
        API.getJob(caseAId),
        API.getJob(caseBId),
      ]);
      if (!jobA.run_dir || !jobB.run_dir) {
        status.textContent = 'One of the selected jobs has no analysis directory.';
        return;
      }

      const payload = {
        case_a: jobA.run_dir,
        case_b: jobB.run_dir,
        label_a: data.get('label_a') || 'Condition A',
        label_b: data.get('label_b') || 'Condition B',
        job_name: data.get('job_name') || '',
        comparison_mode: data.get('comparison_mode'),
        comparison_scope: data.get('comparison_scope'),
        alignment_selection: data.get('alignment_selection'),
        residue_mapping: data.get('residue_mapping'),
        comparison_context: data.get('comparison_context') || '',
        auto_start: form.querySelector('input[name="auto_start"]').checked,
      };

      const job = await API.createCompareJob(payload);
      status.textContent = `Created comparison ${job.id}. Redirecting...`;
      setTimeout(() => { location.hash = `#/jobs/${job.id}`; }, 400);
    } catch (error) {
      status.textContent = `Failed: ${error.message}`;
    }
  });
}

// ========================================================================
// Compare detail page — multi-tab browser for comparison artifacts
// ========================================================================

function deltaBadge(delta, relPct, unit) {
  if (delta === null || delta === undefined || delta === '' || delta === 'NaN' || isNaN(parseFloat(delta))) {
    return '<span style="color: var(--muted); font-size: 12px;">—</span>';
  }
  const d = parseFloat(delta);
  const sign = d > 0 ? '+' : '';
  const color = Math.abs(d) < 1e-6 ? '#6f746d' : (d > 0 ? '#15803d' : '#b91c1c');
  const rel = (relPct !== null && relPct !== '' && !isNaN(parseFloat(relPct))) ? ` (${sign}${parseFloat(relPct).toFixed(1)}%)` : '';
  return `<span style="color: ${color}; font-family: ui-monospace, monospace; font-weight: 700;">${sign}${d.toFixed(3)}${unit ? ' ' + unit : ''}${rel}</span>`;
}

async function renderCompareDetail(jobId, job) {
  let summaryResp, tableResp, plotsResp;
  try {
    [summaryResp, tableResp, plotsResp] = await Promise.all([
      API.getCompareSummary(jobId),
      API.getCompareTable(jobId),
      API.listComparePlots(jobId),
    ]);
  } catch (error) {
    return `<div class="panel"><h2>Comparison not available</h2><p class="muted">${escapeHtml(error.message)}</p><a class="button" href="#/jobs">Back to jobs</a></div>`;
  }

  // Job still running / failed
  if (summaryResp.status === 'no_artifacts') {
    return `
      <div class="panel" style="text-align: center; padding: 40px;">
        <h2>${escapeHtml(job.name || jobId)}</h2>
        <p class="muted">${escapeHtml(summaryResp.message || 'Comparison artifacts not generated yet.')}</p>
        <p style="font-size: 12px; color: var(--muted);">Job status: <strong>${escapeHtml(summaryResp.job_status)}</strong></p>
        <div class="actions" style="justify-content: center;">
          <a class="button secondary" href="#/jobs">Back to jobs</a>
        </div>
      </div>
    `;
  }

  const summary = summaryResp.summary || {};
  const caseA = summary.case_a || {};
  const caseB = summary.case_b || {};
  const comparability = summary.comparability || {};
  const takeaways = summary.takeaways || [];
  const categories = tableResp.categories || {};
  const plots = plotsResp.plots || [];
  const htmlReports = plotsResp.html_reports || [];

  // Comparability badge color
  const compStatus = (comparability.status || '').toLowerCase();
  const compBadge = compStatus === 'comparable'
    ? { bg: '#dbe9e1', color: '#255b4b', label: 'Comparable' }
    : compStatus === 'caution'
      ? { bg: '#f3e8cc', color: '#a36b00', label: 'Caution' }
      : { bg: '#f0ddd5', color: '#9b3a2a', label: comparability.status || 'Unknown' };

  // Tab structure: category → which plots to show
  const tabs = [
    { id: 'overview', label: 'Overview', categories: [] },
    { id: 'quality', label: 'Quality', categories: ['quality'], plots: ['quality_interface_comparison'] },
    { id: 'interface', label: 'Interface', categories: ['interface', 'interaction_family'], plots: ['interaction_family_comparison'] },
    { id: 'flexibility', label: 'Flexibility', categories: ['rmsf'], plots: ['flexibility_comparison'] },
    { id: 'hotspots', label: 'Hotspots (RRCS)', categories: ['rrcs'], plots: ['rrcs_region_comparison'] },
    { id: 'identity', label: 'Identity', categories: ['identity'] },
  ];

  // Filter tabs that have data
  const tabsWithData = tabs.filter(t => {
    if (t.id === 'overview') return true;
    if (t.categories && t.categories.some(c => categories[c] && categories[c].length)) return true;
    if (t.plots && t.plots.some(p => plots.some(pl => pl.name.startsWith(p)))) return true;
    return false;
  });

  return `
    <div style="max-width: 1200px; margin: 0 auto;">
      <!-- Hero -->
      <div style="background: linear-gradient(135deg, #1c473c 0%, #2d6f5c 100%); border-radius: 24px; padding: 28px 32px; color: #f7f5ef; box-shadow: 0 18px 44px rgba(28, 71, 60, 0.18); position: relative; overflow: hidden;">
        <div style="position: absolute; top: -40px; right: -40px; width: 200px; height: 200px; border-radius: 50%; background: radial-gradient(circle, rgba(255,209,102,0.15) 0%, transparent 70%); pointer-events: none;"></div>
        <div style="position: relative;">
          <div style="display: flex; align-items: center; gap: 10px; margin-bottom: 6px;">
            <span style="padding: 3px 10px; border-radius: 999px; background: rgba(255, 209, 102, 0.2); color: #ffd166; font-size: 10px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase;">COMPARISON</span>
            <span style="padding: 3px 10px; border-radius: 999px; background: ${compBadge.bg}; color: ${compBadge.color}; font-size: 10px; font-weight: 800; letter-spacing: 0.08em; text-transform: uppercase;">${escapeHtml(compBadge.label)}</span>
            ${summary.comparison_mode ? `<span style="opacity: 0.7; font-size: 11px;">mode: ${escapeHtml(summary.comparison_mode)}</span>` : ''}
          </div>
          <h1 style="margin: 4px 0; font-size: 26px; font-weight: 800;">${escapeHtml(job.name || jobId)}</h1>
          ${summary.comparison_context ? `<p style="margin: 4px 0; opacity: 0.85; font-size: 13px;">${escapeHtml(summary.comparison_context)}</p>` : ''}

          <div style="display: grid; grid-template-columns: 1fr auto 1fr; gap: 16px; align-items: center; margin-top: 20px;">
            <div style="background: rgba(247, 245, 239, 0.12); backdrop-filter: blur(10px); border-radius: 14px; padding: 14px 18px;">
              <div style="font-size: 10px; letter-spacing: 0.12em; text-transform: uppercase; opacity: 0.85; font-weight: 700;">Case A · reference</div>
              <div style="font-size: 18px; font-weight: 800; margin-top: 6px;">${escapeHtml(caseA.label || 'A')}</div>
              <div style="font-size: 11px; opacity: 0.75; margin-top: 4px; font-family: ui-monospace, monospace; word-break: break-all;">${escapeHtml(caseA.case_root || '')}</div>
            </div>
            <div style="font-size: 22px; opacity: 0.6;">vs</div>
            <div style="background: rgba(247, 245, 239, 0.12); backdrop-filter: blur(10px); border-radius: 14px; padding: 14px 18px;">
              <div style="font-size: 10px; letter-spacing: 0.12em; text-transform: uppercase; opacity: 0.85; font-weight: 700;">Case B · comparison</div>
              <div style="font-size: 18px; font-weight: 800; margin-top: 6px;">${escapeHtml(caseB.label || 'B')}</div>
              <div style="font-size: 11px; opacity: 0.75; margin-top: 4px; font-family: ui-monospace, monospace; word-break: break-all;">${escapeHtml(caseB.case_root || '')}</div>
            </div>
          </div>
        </div>
      </div>

      <!-- Takeaways -->
      ${takeaways.length ? `
        <div class="panel" style="margin-top: 20px; padding: 18px 22px;">
          <span style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: var(--accent-strong);">Key takeaways</span>
          <ul style="margin: 8px 0 0; padding-left: 22px; line-height: 1.7; font-size: 14px;">
            ${takeaways.map(t => `<li>${escapeHtml(t)}</li>`).join('')}
          </ul>
          ${comparability.reasons && comparability.reasons.length ? `
            <div style="margin-top: 12px; padding-top: 10px; border-top: 1px solid var(--line); font-size: 12px; color: var(--muted);">
              <strong>Comparability:</strong> ${comparability.reasons.map(escapeHtml).join(' · ')}
            </div>
          ` : ''}
        </div>
      ` : ''}

      <!-- Tabs -->
      <div style="margin-top: 20px; border-bottom: 1px solid var(--line); display: flex; gap: 4px; overflow-x: auto;">
        ${tabsWithData.map((t, i) => `
          <button class="compare-tab" data-tab="${t.id}" style="padding: 10px 16px; background: ${i === 0 ? 'white' : 'transparent'}; border: 0; border-bottom: 2px solid ${i === 0 ? 'var(--accent)' : 'transparent'}; color: ${i === 0 ? 'var(--accent-strong)' : 'var(--muted)'}; font-weight: 700; font-size: 13px; cursor: pointer; white-space: nowrap; transition: all 150ms ease;">
            ${escapeHtml(t.label)}
          </button>
        `).join('')}
      </div>

      <!-- Tab content -->
      <div id="compare-tab-content" style="margin-top: 16px;">
        ${renderCompareTab(tabsWithData[0], summary, categories, plots, caseA, caseB, jobId)}
      </div>

      <!-- Footer with full report links -->
      ${htmlReports.length ? `
        <div style="margin-top: 24px; padding: 16px; background: rgba(245, 241, 232, 0.5); border-radius: 14px; display: flex; align-items: center; justify-content: space-between; gap: 12px;">
          <div>
            <strong style="font-size: 13px;">Want the full HTML report?</strong>
            <p class="muted" style="margin: 4px 0 0; font-size: 12px;">Open the standalone report with all plots and detail tables.</p>
          </div>
          <div style="display: flex; gap: 8px;">
            ${htmlReports.map(r => `<a class="button" href="${escapeHtml(r.url)}" target="_blank">${escapeHtml(r.label)}</a>`).join('')}
          </div>
        </div>
      ` : ''}

      <!-- Stash data for tab switcher -->
      <script id="compare-data" type="application/json">${JSON.stringify({
        summary,
        categories,
        plots,
        caseA,
        caseB,
        jobId,
        tabsWithData,
      }).replace(/</g, '\\u003c')}</script>
    </div>
  `;
}

function renderCompareTab(tab, summary, categories, plots, caseA, caseB, jobId) {
  if (!tab) return '<div class="panel"><p class="muted">No data.</p></div>';

  if (tab.id === 'overview') {
    return renderCompareOverview(summary, categories, plots, caseA, caseB);
  }

  // Generic tab: show category tables + relevant plots
  const tabPlots = (tab.plots || []).flatMap(prefix =>
    plots.filter(p => p.name.startsWith(prefix))
  );

  const tabCategories = (tab.categories || []);
  const tables = tabCategories.map(cat => ({ cat, rows: categories[cat] || [] })).filter(t => t.rows.length);

  return `
    ${tabPlots.length ? `
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 16px; margin-bottom: 16px;">
        ${tabPlots.map(p => `
          <div class="panel" style="padding: 14px; text-align: center;">
            <div style="font-size: 11px; font-weight: 800; color: var(--muted); letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 8px;">${escapeHtml(p.label)}</div>
            <img src="${escapeHtml(p.url)}" alt="${escapeHtml(p.label)}" style="max-width: 100%; max-height: 400px; height: auto; border-radius: 8px;">
          </div>
        `).join('')}
      </div>
    ` : ''}

    ${tables.map(t => renderCategoryTable(t.cat, t.rows, caseA.label, caseB.label)).join('') ||
      '<div class="panel"><p class="muted">No detailed metrics available in this section.</p></div>'}
  `;
}

function renderCompareOverview(summary, categories, plots, caseA, caseB) {
  // Show all 4 main plots in a grid + a "metric overview" table snapshot
  const overviewPlots = plots.filter(p =>
    p.name.includes('quality_interface') ||
    p.name.includes('flexibility') ||
    p.name.includes('rrcs_region') ||
    p.name.includes('interaction_family')
  );
  const metricCount = Object.values(categories).reduce((acc, rows) => acc + rows.length, 0);

  return `
    <!-- Quick stats -->
    <div class="grid cards" style="margin-bottom: 16px; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));">
      <div class="card" style="padding: 14px;">
        <span>Metrics</span><strong style="font-size: 26px;">${metricCount}</strong>
      </div>
      <div class="card" style="padding: 14px;">
        <span>Categories</span><strong style="font-size: 26px;">${Object.keys(categories).length}</strong>
      </div>
      <div class="card" style="padding: 14px;">
        <span>Plots</span><strong style="font-size: 26px;">${plots.length}</strong>
      </div>
    </div>

    <!-- Plot grid -->
    ${overviewPlots.length ? `
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); gap: 14px;">
        ${overviewPlots.map(p => `
          <div class="panel" style="padding: 12px; text-align: center;">
            <div style="font-size: 11px; font-weight: 800; color: var(--muted); letter-spacing: 0.1em; text-transform: uppercase; margin-bottom: 8px;">${escapeHtml(p.label)}</div>
            <img src="${escapeHtml(p.url)}" alt="${escapeHtml(p.label)}" style="max-width: 100%; max-height: 280px; height: auto; border-radius: 6px;">
          </div>
        `).join('')}
      </div>
    ` : '<div class="panel"><p class="muted">No plots generated yet.</p></div>'}
  `;
}

function renderCategoryTable(category, rows, labelA, labelB) {
  if (!rows.length) return '';
  const titleMap = {
    quality: 'Quality & convergence',
    interface: 'Interface (BSA, contacts)',
    rmsf: 'Flexibility (RMSF by region)',
    rrcs: 'RRCS hotspots (by region)',
    interaction_family: 'Interaction families',
    identity: 'System identity',
  };
  return `
    <div class="panel" style="margin-bottom: 14px;">
      <h3 style="margin: 0 0 12px;">${escapeHtml(titleMap[category] || category)}</h3>
      <div style="overflow-x: auto;">
        <table style="width: 100%; font-size: 13px;">
          <thead>
            <tr>
              <th>Metric</th>
              <th style="text-align: right;">${escapeHtml(labelA || 'A')}</th>
              <th style="text-align: right;">${escapeHtml(labelB || 'B')}</th>
              <th style="text-align: right;">Δ (B − A)</th>
            </tr>
          </thead>
          <tbody>
            ${rows.map(r => `
              <tr>
                <td><strong>${escapeHtml(r.metric)}</strong>${r.unit ? `<span class="muted" style="font-size: 11px; margin-left: 6px;">${escapeHtml(r.unit)}</span>` : ''}</td>
                <td style="text-align: right; font-family: ui-monospace, monospace;">${escapeHtml(r.case_a || '—')}</td>
                <td style="text-align: right; font-family: ui-monospace, monospace;">${escapeHtml(r.case_b || '—')}</td>
                <td style="text-align: right;">${deltaBadge(r.delta, r.relative_change_percent, r.unit)}</td>
              </tr>
            `).join('')}
          </tbody>
        </table>
      </div>
    </div>
  `;
}

function mountCompareDetail(jobId) {
  // Tab switching
  const dataEl = document.getElementById('compare-data');
  if (!dataEl) return;
  let bundle;
  try { bundle = JSON.parse(dataEl.textContent); } catch (e) { return; }
  const { summary, categories, plots, caseA, caseB, tabsWithData } = bundle;
  const content = document.getElementById('compare-tab-content');

  document.querySelectorAll('.compare-tab').forEach(btn => {
    btn.addEventListener('click', () => {
      const tabId = btn.getAttribute('data-tab');
      const tab = tabsWithData.find(t => t.id === tabId);
      if (!tab || !content) return;
      // Update visual state
      document.querySelectorAll('.compare-tab').forEach(b => {
        b.style.background = 'transparent';
        b.style.borderBottom = '2px solid transparent';
        b.style.color = 'var(--muted)';
      });
      btn.style.background = 'white';
      btn.style.borderBottom = '2px solid var(--accent)';
      btn.style.color = 'var(--accent-strong)';
      // Render content
      content.innerHTML = renderCompareTab(tab, summary, categories, plots, caseA, caseB, jobId);
    });
  });
}

async function renderNewAnalysis() {
  const reuseJobId = window.App.queryParams?.get('reuse');
  let reuseBanner = '';
  let reuseInfo = null;

  if (reuseJobId) {
    try {
      const trajectories = await API.listTrajectories();
      reuseInfo = trajectories.find(t => t.job_id === reuseJobId);
      if (reuseInfo) {
        reuseBanner = `
          <div style="margin-bottom: 16px; padding: 14px 16px; border-radius: 14px; background: linear-gradient(90deg, rgba(20, 184, 166, 0.1), rgba(14, 165, 233, 0.06)); border: 1px solid rgba(20, 184, 166, 0.24); display: flex; justify-content: space-between; align-items: center; gap: 12px;">
            <div>
              <div style="font-size: 11px; font-weight: 800; color: var(--accent-strong); text-transform: uppercase; letter-spacing: 0.06em; margin-bottom: 4px;">Reusing trajectory</div>
              <div style="font-size: 14px; font-weight: 700; color: var(--ink);">
                <span style="font-family: ui-monospace, monospace;">${escapeHtml(reuseInfo.job_name || reuseInfo.job_id)}</span>
                <span style="color: var(--muted); font-weight: 500; margin-left: 8px;">${reuseInfo.quality?.duration_ns || '-'} ns · ${reuseInfo.chain_info?.num_chains || 0} chains · stride ${reuseInfo.stride}</span>
              </div>
              <div style="font-size: 12px; color: var(--muted); margin-top: 4px;">PBC correction will be skipped — analysis modules run directly on the prepared trajectory.</div>
            </div>
            <a href="#/new" style="padding: 6px 12px; border-radius: 10px; background: white; border: 1px solid var(--line); color: var(--muted); text-decoration: none; font-size: 12px; font-weight: 700; white-space: nowrap;">Cancel reuse</a>
          </div>
        `;
      }
    } catch (e) {
      console.warn('Failed to load trajectory info:', e);
    }
  }

  return `
    <form id="new-job-form" class="panel analysis-form">
      ${reuseBanner}
      ${reuseJobId ? `<input type="hidden" name="reuse_from_job_id" value="${escapeHtml(reuseJobId)}">` : ''}
      <div class="form-header">
        <div>
          <h2>${reuseInfo ? 'Reuse Trajectory' : 'Create Analysis Job'}</h2>
          <p class="muted">${reuseInfo ? 'Skip preprocessing and run analysis modules directly on a previously prepared trajectory.' : 'Upload raw or prepared trajectory inputs. ImmunoScope will validate the inputs and execute the selected modules through the CLI workflow.'}</p>
        </div>
        <label class="switch-row">
          <input name="auto_start" type="checkbox" checked value="true">
          Start immediately
        </label>
      </div>

      <div class="grid two-col">
        <label>
          Job name
          <input name="job_name" placeholder="1mi5_standard_analysis">
        </label>
        <label>
          Profile
          <select name="profile">
            <option value="standard">Standard single-trajectory analysis</option>
            <option value="custom">Custom workflow</option>
          </select>
        </label>
      </div>

      <div class="grid two-col">
        <label>
          Input mode
          <select name="input_mode" id="input-mode">
            <option value="raw">Raw trajectory: topology + trajectory (structure optional)</option>
            <option value="prepared">Prepared trajectory: prepared structure + topology + processed trajectory</option>
          </select>
        </label>
        <label>
          Module preset
          <select id="module-preset">
            <option value="standard">Standard full report</option>
            <option value="contact">Contact/RRCS report</option>
            <option value="geometry">Geometry (angles + dihedrals)</option>
            <option value="fast">Fast intake report</option>
            <option value="custom">Custom modules</option>
          </select>
        </label>
      </div>

      <section class="input-block" data-mode="raw">
        <h3>Raw Inputs</h3>
        <div class="grid three-col">
          ${renderFileField('structure', 'Structure PDB (optional)', '.pdb,.gro', false)}
          ${renderFileField('topology', 'Topology TPR', '.tpr')}
          ${renderFileField('trajectory', 'Trajectory XTC/TRR', '.xtc,.trr')}
        </div>
        <p class="muted" style="margin-top: 8px; font-size: 12px;">💡 Structure is optional - if not provided, it will be extracted from the TPR file automatically.</p>
      </section>

      <section class="input-block hidden" data-mode="prepared">
        <h3>Prepared Inputs</h3>
        <div class="grid three-col">
          ${renderFileField('prepared_structure', 'Prepared structure', '.pdb,.gro')}
          ${renderFileField('topology_prepared', 'Topology TPR', '.tpr')}
          ${renderFileField('processed_trajectory', 'Processed trajectory', '.xtc,.trr')}
        </div>
      </section>

      <div class="grid three-col">
        <label>
          Preprocess method
          <select name="preprocess_method">
            <option value="2step">2step</option>
            <option value="3step">3step</option>
          </select>
        </label>
        <label>Stride<input name="stride" type="number" min="1" value="5"></label>
        <label>
          Landscape reducer
          <select name="landscape_reducer">
            <option value="pca">PCA</option>
            <option value="umap">UMAP</option>
            <option value="tica">TICA</option>
          </select>
        </label>
      </div>

      <section class="module-picker">
        <div class="module-picker-head">
          <div>
            <h3>Analysis Modules</h3>
            <p class="muted">Select the analysis steps to include in this job.</p>
          </div>
          <span id="module-count" class="module-count"></span>
        </div>
        <div class="module-grid">
          ${ANALYSIS_MODULES.map(module => renderModuleOption(module)).join('')}
        </div>
        <input type="hidden" name="modules" id="modules-input" value="${MODULE_PRESETS.standard}">
      </section>
      <label>
        Notes
        <textarea name="notes" placeholder="Optional context for this job"></textarea>
      </label>
      <button type="submit">Create job</button>
      <p id="new-job-status" class="muted"></p>
    </form>
  `;
}

function renderModuleOption(module) {
  return `
    <label class="module-option">
      <input type="checkbox" value="${escapeHtml(module.id)}">
      <span class="module-option-body">
        <span class="module-option-title">${escapeHtml(module.label)}</span>
        <span class="module-option-description">${escapeHtml(module.description)}</span>
      </span>
    </label>
  `;
}

function renderFileField(name, label, accept, required = true) {
  const requiredAttr = required ? 'required' : '';
  return `
    <label class="file-field">
      <span class="field-label">${escapeHtml(label)}</span>
      <input name="${escapeHtml(name)}" type="file" accept="${escapeHtml(accept)}" ${requiredAttr}>
      <span class="file-picker">
        <span class="file-picker-title">Select file</span>
        <span class="file-picker-name" data-empty="No file selected">No file selected</span>
      </span>
    </label>
  `;
}

function mountNewAnalysis() {
  const form = document.getElementById('new-job-form');
  const status = document.getElementById('new-job-status');
  const mode = document.getElementById('input-mode');
  const preset = document.getElementById('module-preset');
  const modules = document.getElementById('modules-input');
  const moduleCount = document.getElementById('module-count');
  const moduleChecks = Array.from(form.querySelectorAll('.module-option input[type="checkbox"]'));

  // Check if in reuse mode - hide file upload sections
  const isReuseMode = !!form.querySelector('input[name="reuse_from_job_id"]');
  if (isReuseMode) {
    document.querySelectorAll('.input-block').forEach(block => {
      block.style.display = 'none';
      block.querySelectorAll('input[type="file"]').forEach(input => input.removeAttribute('required'));
    });
    // Also hide the input mode selector
    const modeContainer = mode?.closest('label');
    if (modeContainer) modeContainer.style.display = 'none';
  }

  function updateMode() {
    if (isReuseMode) return;  // Skip mode toggling in reuse mode
    document.querySelectorAll('.input-block').forEach(block => {
      const isHidden = block.dataset.mode !== mode.value;
      block.classList.toggle('hidden', isHidden);

      // Toggle required attribute on file inputs based on visibility
      block.querySelectorAll('input[type="file"]').forEach(input => {
        if (isHidden) {
          input.removeAttribute('required');
        } else {
          // Only set required for non-optional fields
          const fieldName = input.getAttribute('name');
          if (fieldName !== 'structure') {  // structure is optional in raw mode
            input.setAttribute('required', 'required');
          }
        }
      });
    });
  }
  function selectedModules() {
    return moduleChecks.filter(input => input.checked).map(input => input.value);
  }
  function applyModules(value) {
    const selected = new Set(String(value || '').split(',').map(item => item.trim()).filter(Boolean));
    moduleChecks.forEach(input => {
      input.checked = selected.has(input.value);
    });
    syncModules();
  }
  function syncModules() {
    const selected = selectedModules();
    modules.value = selected.join(',');
    moduleCount.textContent = `${selected.length} selected`;
  }
  function syncPresetFromSelection() {
    const current = selectedModules().join(',');
    const matched = Object.entries(MODULE_PRESETS).find(([, value]) => value === current);
    preset.value = matched ? matched[0] : 'custom';
  }
  mode.addEventListener('change', updateMode);
  preset.addEventListener('change', () => {
    if (preset.value !== 'custom') applyModules(MODULE_PRESETS[preset.value]);
  });
  moduleChecks.forEach(input => {
    input.addEventListener('change', () => {
      syncModules();
      syncPresetFromSelection();
    });
  });
  applyModules(MODULE_PRESETS.standard);
  updateMode();  // Initialize mode on page load
  form.querySelectorAll('.file-field input[type="file"]').forEach(input => {
    input.addEventListener('change', () => {
      const name = input.files && input.files.length ? input.files[0].name : 'No file selected';
      input.closest('.file-field').querySelector('.file-picker-name').textContent = name;
    });

    // Add drag-and-drop support
    const fileField = input.closest('.file-field');
    const filePicker = fileField.querySelector('.file-picker');

    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
      filePicker.addEventListener(eventName, preventDefaults, false);
    });

    function preventDefaults(e) {
      e.preventDefault();
      e.stopPropagation();
    }

    ['dragenter', 'dragover'].forEach(eventName => {
      filePicker.addEventListener(eventName, () => {
        filePicker.classList.add('drag-over');
      }, false);
    });

    ['dragleave', 'drop'].forEach(eventName => {
      filePicker.addEventListener(eventName, () => {
        filePicker.classList.remove('drag-over');
      }, false);
    });

    filePicker.addEventListener('drop', (e) => {
      const dt = e.dataTransfer;
      const files = dt.files;
      if (files.length > 0) {
        input.files = files;
        const name = files[0].name;
        filePicker.querySelector('.file-picker-name').textContent = name;
      }
    }, false);
  });

  form.addEventListener('submit', async event => {
    event.preventDefault();
    status.textContent = 'Creating job...';
    const data = new FormData(form);
    const autoStart = form.querySelector('input[name="auto_start"]').checked;
    data.set('auto_start', autoStart ? 'true' : 'false');
    if (mode.value === 'prepared') {
      const preparedTopology = form.querySelector('input[name="topology_prepared"]').files[0];
      if (preparedTopology) data.set('topology', preparedTopology);
      data.delete('structure');
      data.delete('trajectory');
    } else {
      data.delete('prepared_structure');
      data.delete('processed_trajectory');
      data.delete('topology_prepared');
    }
    try {
      const job = await API.createJob(data);
      status.textContent = `Created job ${job.id}.`;
      location.hash = `#/jobs/${job.id}`;
    } catch (error) {
      status.textContent = `Failed: ${error.message}`;
    }
  });
}

async function renderJobs() {
  const jobs = await API.listJobs().catch(() => []);
  if (!jobs.length) {
    return `
      <div class="panel">
        <h2>No jobs yet</h2>
        <p class="muted">Create an analysis job to populate this table.</p>
        <a class="button" href="#/new">Create job</a>
      </div>
    `;
  }
  const rows = jobs.map(job => `
    <tr>
      <td><a href="#/jobs/${escapeHtml(job.id)}"><strong>${escapeHtml(job.name)}</strong></a><br><span class="muted">${escapeHtml(job.id)}</span></td>
      <td>${statusBadge(job.status)}<div class="progress"><span style="width:${Number(job.progress || 0)}%"></span></div></td>
      <td>${escapeHtml(job.config.profile)}<br><span class="muted">${escapeHtml(job.config.input_mode)}</span></td>
      <td>${escapeHtml((job.config.modules || []).join(', '))}</td>
      <td>${formatTime(job.created_at)}</td>
    </tr>
  `).join('');
  return `
    <div class="panel">
      <div class="section-header">
        <h2>Jobs</h2>
        <a class="button secondary" href="#/new">New job</a>
      </div>
      <table>
        <thead><tr><th>Job</th><th>Status</th><th>Profile</th><th>Modules</th><th>Created</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </div>
  `;
}

async function renderJobDetail(id) {
  const job = await API.getJob(id).catch(error => ({ error: error.message }));
  if (job.error && !job.id) {
    return `<div class="panel"><h2>Job not available</h2><p class="muted">${escapeHtml(job.error)}</p></div>`;
  }
  const summary = job.run_summary || {};
  const stages = ((summary.manifest && []) || []);
  const index = job.result_index_payload || {};
  const moduleRows = (index.modules || []).map(module => `
    <tr>
      <td>${escapeHtml(module.name)}</td>
      <td>${statusBadge(module.status)}</td>
      <td>${escapeHtml(Object.entries(module.metrics || {}).map(([k, v]) => `${k}: ${v}`).join(', '))}</td>
    </tr>
  `).join('');
  return `
    <div class="job-layout">
      <section class="panel">
        <div class="section-header">
          <div>
            <h2>${escapeHtml(job.name)}</h2>
            <p class="muted">${escapeHtml(job.id)} | ${escapeHtml(job.config.input_mode)} | ${escapeHtml(job.config.modules.join(', '))}</p>
          </div>
          <div class="actions">
            ${job.status === 'queued' || job.status === 'failed' ? `<button id="start-job">Start</button>` : ''}
            ${job.report_html ? `<a class="button secondary" href="/api/jobs/${escapeHtml(job.id)}/report" target="_blank">Open report</a>` : ''}
          </div>
        </div>
        <div class="grid cards compact-cards">
          <div class="card"><span>Status</span><strong>${escapeHtml(job.status)}</strong></div>
          <div class="card"><span>Progress</span><strong>${Number(job.progress || 0)}%</strong></div>
          <div class="card"><span>Started</span><strong>${escapeHtml(formatTime(job.started_at))}</strong></div>
          <div class="card"><span>Finished</span><strong>${escapeHtml(formatTime(job.finished_at))}</strong></div>
        </div>
        ${job.error ? `<div class="alert">${escapeHtml(job.error)}</div>` : ''}
      </section>

      <section class="panel">
        <h3>Inputs</h3>
        <div class="kv">
          ${Object.entries(job.inputs || {}).map(([key, value]) => `<span>${escapeHtml(key)}</span><code>${escapeHtml(value || '-')}</code>`).join('')}
        </div>
      </section>

      <section class="panel">
        <h3>Module Results</h3>
        ${moduleRows ? `<table><thead><tr><th>Module</th><th>Status</th><th>Metrics</th></tr></thead><tbody>${moduleRows}</tbody></table>` : '<p class="muted">No module result index is available yet.</p>'}
      </section>

      <section class="panel assistant-panel">
        <h3>Query-Agent <span style="font-size:11px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:0.06em;margin-left:6px;">read-only · this job</span></h3>
        <div id="assistant-skill-strip" class="assistant-skill-strip">
          <span class="muted">Loading available skills...</span>
        </div>
        <div class="quick-actions">
          ${QUICK_QUESTIONS.map(q => `<button class="button secondary quick-question" data-question="${escapeHtml(q)}">${escapeHtml(q)}</button>`).join('')}
        </div>
        <div id="assistant-messages" class="assistant-messages">
          <div class="assistant-message assistant">Read-only Q&amp;A about <strong>this job's</strong> results — it queries the analysis you already ran. For mutation <em>design</em> (multi-agent readers + critic + saved picks), use <a href="#/copilot/new">ImmunoScope</a>.</div>
        </div>
        <form id="assistant-form" class="assistant-form">
          <input id="assistant-input" placeholder="Ask about stability, interface, hotspots, or design sites">
          <button type="submit">Send</button>
        </form>
      </section>

      <section class="panel">
        <h3>Log</h3>
        <pre class="log-view">${escapeHtml(job.log || 'No log yet.')}</pre>
      </section>
    </div>
  `;
}

function mountJobDetail(id) {
  const startButton = document.getElementById('start-job');
  if (startButton) {
    startButton.addEventListener('click', async () => {
      startButton.disabled = true;
      await API.startJob(id).catch(() => null);
      await App.navigate();
    });
  }
  const form = document.getElementById('assistant-form');
  if (!form) return;
  const input = document.getElementById('assistant-input');
  const messages = document.getElementById('assistant-messages');
  const skillStrip = document.getElementById('assistant-skill-strip');
  API.listJobSkills(id).then(payload => {
    const skills = payload.skills || [];
    if (!skills.length) {
      skillStrip.innerHTML = '<span class="muted">No result-backed skills are available yet.</span>';
      return;
    }
    skillStrip.innerHTML = skills.map(skill => `
      <button class="skill-pill" data-skill="${escapeHtml(skill.name)}" title="${escapeHtml(skill.description || '')}">
        ${escapeHtml(skill.name)}
      </button>
    `).join('');
    skillStrip.querySelectorAll('.skill-pill').forEach(button => {
      button.addEventListener('click', async () => {
        await send(button.dataset.skill || '', button.dataset.skill || '');
      });
    });
  }).catch(error => {
    skillStrip.innerHTML = `<span class="muted">Skill list unavailable: ${escapeHtml(error.message)}</span>`;
  });
  async function send(question, skill = null) {
    if (!question.trim()) return;
    messages.insertAdjacentHTML('beforeend', `<div class="assistant-message user">${escapeHtml(question)}</div>`);
    const loading = document.createElement('div');
    loading.className = 'assistant-message assistant';
    loading.textContent = 'Thinking...';
    messages.appendChild(loading);
    input.value = '';
    try {
      const result = skill ? await API.invokeJobSkill(id, { skill, question }) : await API.askJob(id, { question });
      const evidenceCards = ((result.agent_result_schema || {}).evidence_cards || []).slice(0, 3);
      loading.innerHTML = `
        <strong>${escapeHtml(result.skill || 'assistant')}</strong>
        <p>${escapeHtml(result.answer || 'No answer available.')}</p>
        ${(result.evidence || []).slice(0, 4).map(item => `<div class="evidence-line">${escapeHtml(item)}</div>`).join('')}
        ${evidenceCards.length ? `<div class="evidence-cards">${evidenceCards.map(card => `
          <div class="mini-card">
            <span>${escapeHtml(card.module)}</span>
            <strong>${escapeHtml(card.summary)}</strong>
          </div>
        `).join('')}</div>` : ''}
      `;
    } catch (error) {
      loading.textContent = `Assistant failed: ${error.message}`;
      loading.classList.add('error');
    }
    messages.scrollTop = messages.scrollHeight;
  }
  form.addEventListener('submit', event => {
    event.preventDefault();
    send(input.value);
  });
  document.querySelectorAll('.quick-question').forEach(button => {
    button.addEventListener('click', () => send(button.dataset.question || ''));
  });
}

// Agent page render
function renderAgent() {
  return `
    <div class="agent-container">
      <div class="agent-header">
        <div class="agent-status">
          <span class="status-dot" id="statusDot"></span>
          <span id="statusText">Connecting...</span>
        </div>
      </div>
      <div class="error-banner" id="errorBanner"></div>
      <div class="chat-container" id="chatContainer">
        <div class="welcome-message">
          <h3>Query-Agent</h3>
          <p>Ask questions across your MD analyses — trajectory quality, the interface, hotspots, or a specific residue. For mutation <em>design</em> (multi-agent readers + critic), use <strong>ImmunoScope</strong>.</p>
          <div class="example-questions">
            <button class="example-btn" onclick="window.App.agentClient.sendMessage('Analyze output/5c0a_run2_full_analysis')">Analyze a trajectory</button>
            <button class="example-btn" onclick="window.App.agentClient.sendMessage('What are the key hotspots?')">Find hotspots</button>
            <button class="example-btn" onclick="location.hash='#/copilot'">Design mutations →</button>
          </div>
        </div>
      </div>
      <div class="input-container">
        <textarea id="messageInput" placeholder="Type your message..." rows="1"></textarea>
        <button id="sendButton" onclick="window.App.agentClient.sendMessage(document.getElementById('messageInput').value)">Send</button>
      </div>
    </div>
  `;
}

// Trajectories page - asset library of processed trajectories
async function renderTrajectories() {
  const trajectories = await API.listTrajectories().catch(() => []);

  if (!trajectories.length) {
    return `
      <div class="panel">
        <h2>Trajectory Library</h2>
        <p class="muted">No processed trajectories yet. Trajectories appear here once preprocessing (PBC correction) completes.</p>
        <div class="actions">
          <a class="button" href="#/new">Upload new trajectory</a>
        </div>
      </div>
    `;
  }

  const formatBytes = (bytes) => {
    if (!bytes) return '-';
    if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(0)} KB`;
    if (bytes < 1024 * 1024 * 1024) return `${(bytes / 1024 / 1024).toFixed(1)} MB`;
    return `${(bytes / 1024 / 1024 / 1024).toFixed(2)} GB`;
  };

  const convergenceBadge = (quality) => {
    const map = {
      excellent: { color: '#15803d', bg: '#dcfce7', label: 'Excellent' },
      good: { color: '#0f766e', bg: '#ccfbf1', label: 'Good' },
      moderate: { color: '#a16207', bg: '#fef3c7', label: 'Moderate' },
      unstable: { color: '#b91c1c', bg: '#fee2e2', label: 'Unstable' },
      unknown: { color: '#64748b', bg: '#f1f5f9', label: 'Unknown' },
    };
    const cfg = map[quality] || map.unknown;
    return `<span style="padding: 3px 8px; border-radius: 999px; background: ${cfg.bg}; color: ${cfg.color}; font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em;">${cfg.label}</span>`;
  };

  const cards = trajectories.map(traj => {
    const chains = Object.keys(traj.chain_info?.chain_mapping || {});
    const q = traj.quality || {};

    return `
      <div class="card" style="display: flex; flex-direction: column; gap: 0; padding: 0; overflow: hidden;">
        <div style="padding: 18px 18px 14px; background: linear-gradient(180deg, rgba(236, 254, 255, 0.4), transparent);">
          <div style="display: flex; justify-content: space-between; align-items: flex-start; gap: 12px;">
            <div style="flex: 1; min-width: 0;">
              <span style="font-size: 10px; font-weight: 800; color: var(--muted); text-transform: uppercase; letter-spacing: 0.1em;">Trajectory</span>
              <strong style="display: block; margin-top: 4px; font-size: 16px; font-family: ui-monospace, monospace; word-break: break-all;">${escapeHtml(traj.job_name || traj.job_id)}</strong>
            </div>
            ${convergenceBadge(q.convergence_quality)}
          </div>
        </div>

        <div style="padding: 14px 18px; display: grid; gap: 10px; font-size: 13px;">
          <div style="display: flex; justify-content: space-between; align-items: center;">
            <span class="muted" style="font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em;">Simulation</span>
            <span style="font-weight: 700; font-family: ui-monospace, monospace;">${q.duration_ns || '-'} ns · ${q.n_frames || '-'} frames</span>
          </div>

          <div style="display: flex; justify-content: space-between; align-items: center;">
            <span class="muted" style="font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em;">Chains</span>
            <span style="font-weight: 700;">${traj.chain_info?.num_chains || 0} chains <span style="color: var(--muted); font-family: ui-monospace, monospace; font-size: 12px;">(${chains.sort().join(', ') || '-'})</span></span>
          </div>

          <div style="display: flex; justify-content: space-between; align-items: center;">
            <span class="muted" style="font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em;">RMSD</span>
            <span style="font-weight: 700; font-family: ui-monospace, monospace;">${q.tail90_mean_rmsd_nm || '-'} ± ${q.tail90_std_rmsd_nm || '-'} nm</span>
          </div>

          <div style="display: flex; justify-content: space-between; align-items: center;">
            <span class="muted" style="font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em;">Size</span>
            <span style="font-weight: 700; font-family: ui-monospace, monospace;">${formatBytes(traj.trajectory_size_bytes)}</span>
          </div>

          <div style="display: flex; justify-content: space-between; align-items: center;">
            <span class="muted" style="font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em;">Preprocess</span>
            <span style="font-weight: 700;">${escapeHtml(traj.preprocess_method)} · stride ${traj.stride}</span>
          </div>
        </div>

        <div style="padding: 12px 18px; border-top: 1px solid var(--line); background: rgba(248, 250, 252, 0.6); font-size: 11px; color: var(--muted);">
          ${formatTime(traj.finished_at || traj.created_at)}
        </div>

        <div style="padding: 12px 18px 18px; display: flex; gap: 8px;">
          ${traj.can_reuse ? `<a class="button" href="#/new?reuse=${escapeHtml(traj.job_id)}" style="flex: 1; justify-content: center;">Reuse</a>` : ''}
          <a class="button secondary" href="#/jobs/${escapeHtml(traj.job_id)}" style="flex: 1; justify-content: center;">Detail</a>
        </div>
      </div>
    `;
  }).join('');

  const reusableCount = trajectories.filter(t => t.can_reuse).length;

  return `
    <div class="panel">
      <div class="form-header">
        <div>
          <h2>Trajectory Library</h2>
          <p class="muted">${trajectories.length} processed trajector${trajectories.length === 1 ? 'y' : 'ies'} · ${reusableCount} ready to reuse</p>
        </div>
        <a class="button" href="#/new">Upload new</a>
      </div>
      <div class="grid cards" style="margin-top: 16px; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 16px;">
        ${cards}
      </div>
    </div>
  `;
}

// Reports page - shows all jobs with generated reports
async function renderReports() {
  const jobs = await API.listJobs().catch(() => []);
  const completedJobs = jobs.filter(job => job.status === 'completed' && job.report_html);

  if (!completedJobs.length) {
    return `
      <div class="panel">
        <h2>Reports</h2>
        <p class="muted">No reports have been generated yet.</p>
        <p>Complete an analysis job to generate a report.</p>
        <div class="actions">
          <a class="button" href="#/new">Create new analysis</a>
          <a class="button secondary" href="#/jobs">View all jobs</a>
        </div>
      </div>
    `;
  }

  const cards = completedJobs.map(job => `
    <div class="card" style="cursor: pointer;" onclick="window.open('/api/jobs/${escapeHtml(job.id)}/report', '_blank')">
      <span>Report</span>
      <strong style="font-size: 18px; margin-top: 12px;">${escapeHtml(job.name || job.id)}</strong>
      <p class="muted" style="margin-top: 8px; font-size: 13px;">
        ${escapeHtml(job.profile || 'standard')} · ${formatTime(job.completed_at || job.updated_at)}
      </p>
      <div class="actions" style="margin-top: 12px;">
        <a class="button" href="/api/jobs/${escapeHtml(job.id)}/report" target="_blank">Open report</a>
        <a class="button secondary" href="#/jobs/${escapeHtml(job.id)}">Job detail</a>
      </div>
    </div>
  `).join('');

  return `
    <div class="panel">
      <div class="form-header">
        <div>
          <h2>Reports</h2>
          <p class="muted">${completedJobs.length} report${completedJobs.length === 1 ? '' : 's'} generated from completed analyses.</p>
        </div>
        <a class="button" href="#/new">New analysis</a>
      </div>
      <div class="grid cards" style="margin-top: 16px;">
        ${cards}
      </div>
    </div>
  `;
}


// ========================================================================
// ImmunoScope (design surface) — Conversational Mutation Design Workspace
// ========================================================================

const DESIGN_GOALS = [
  { id: 'affinity_enhancement', label: 'Affinity Enhancement', description: 'Strengthen TCR-pMHC binding affinity' },
  { id: 'stability', label: 'Stability', description: 'Improve thermodynamic / kinetic stability' },
  { id: 'specificity', label: 'Specificity', description: 'Enhance peptide selectivity' },
  { id: 'cross_reactivity', label: 'Cross-reactivity', description: 'Tune cross-reactivity profile' },
];

// Wire-level task keys — must stay in sync with
// `immunoscope.web.routers.design.SUPPORTED_TASK_KEYS`. Sending anything
// not in this list will be rejected by the router's pydantic validator.
const DESIGN_TASKS = [
  {
    id: 'tcr_affinity',
    label: 'TCR Affinity',
    description: 'Design TCR mutations against a fixed peptide-HLA target.',
  },
  {
    id: 'peptide_presentation',
    label: 'Peptide Presentation',
    description: 'Design peptide mutations to tune HLA presentation stability (anchor-aware).',
  },
];

const DESIGN_TASK_LABEL = Object.fromEntries(DESIGN_TASKS.map(t => [t.id, t.label]));

const PRIORITY_COLORS = {
  high: { bg: '#dcfce7', color: '#15803d' },
  medium: { bg: '#fef3c7', color: '#a16207' },
  low: { bg: '#f1f5f9', color: '#64748b' },
};

// At-a-glance domain label for a recommendation, e.g. "TCRβ · CDR3".
// `chain` is already friendly ("TCRβ") and `region` is the CDR loop; either
// may be missing, so join only what's present.
function chainRegionLabel(d) {
  const parts = [];
  if (d && d.chain) parts.push(d.chain);
  if (d && d.region) parts.push(d.region);
  return parts.join(' · ');
}

// A site whose only "candidates" are its own wild type is a preserve /
// do-not-mutate call, not a mutation — never render it as "R → R".
function isPreserveDraft(d) {
  if (!d) return false;
  if (d.is_preserve) return true;
  const wt = (d.current_aa || '').toUpperCase();
  const muts = (d.suggested_mutations || []).filter(m => m && m.toUpperCase() !== wt);
  return muts.length === 0;
}

// ---------- List page ----------

async function renderDesignCopilot() {
  const [tasks, jobs] = await Promise.all([
    API.listDesignJobs().catch(() => []),
    API.listJobs().catch(() => []),
  ]);
  const completedJobs = jobs.filter(j => j.status === 'completed');

  if (!completedJobs.length) {
    return `
      <div class="panel">
        <h2>ImmunoScope</h2>
        <p class="muted">No completed analyses yet. ImmunoScope works on top of a completed analysis — it needs MD evidence to reason about.</p>
        <a class="button" href="#/new">Run an analysis first</a>
      </div>
    `;
  }

  const taskCards = tasks.length ? tasks.map(t => `
    <div class="card" style="display: flex; flex-direction: column; gap: 10px;">
      <div style="display: flex; justify-content: space-between; align-items: flex-start; gap: 12px;">
        <div style="flex: 1; min-width: 0;">
          <span style="font-size: 10px; font-weight: 800; color: var(--muted); text-transform: uppercase; letter-spacing: 0.1em;">Design conversation</span>
          <strong style="display: block; margin-top: 4px; font-size: 15px;">${escapeHtml(t.name)}</strong>
          <p class="muted" style="margin: 4px 0 0; font-size: 11px;">
            <span style="display: inline-block; padding: 2px 6px; margin-right: 6px; border-radius: 4px; background: var(--accent-soft); color: var(--accent-strong); font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em;">${escapeHtml(DESIGN_TASK_LABEL[t.task] || t.task || 'tcr affinity')}</span>
            From <span style="font-family: ui-monospace, monospace;">${escapeHtml(t.source_job_name)}</span>
          </p>
        </div>
        <span style="padding: 4px 10px; border-radius: 999px; background: rgba(20, 184, 166, 0.1); color: var(--accent-strong); font-size: 10px; font-weight: 800; text-transform: uppercase;">${escapeHtml(t.status)}</span>
      </div>

      ${t.design_goals && t.design_goals.length ? `
        <div style="display: flex; flex-wrap: wrap; gap: 6px;">
          ${t.design_goals.map(g => `<span style="padding: 3px 8px; border-radius: 6px; background: var(--accent-soft); color: var(--accent-strong); font-size: 11px; font-weight: 700;">${escapeHtml(String(g).replace(/_/g, ' '))}</span>`).join('')}
        </div>
      ` : ''}

      <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 8px; padding: 10px 0; border-top: 1px solid var(--line); font-size: 12px;">
        <div><span class="muted">Saved drafts</span><br><strong style="font-size: 16px;">${t.n_drafts || 0}</strong></div>
        <div><span class="muted">Updated</span><br><strong style="font-size: 11px; font-weight: 600;">${formatTime(t.updated_at)}</strong></div>
      </div>

      <div class="actions" style="margin-top: 4px;">
        <a class="button" href="#/copilot/${escapeHtml(t.id)}" style="flex: 1; justify-content: center;">Open Conversation</a>
        <button data-delete="${escapeHtml(t.id)}" class="button secondary" style="padding: 11px 14px;">Delete</button>
      </div>
    </div>
  `).join('') : '<div class="panel" style="text-align: center;"><p class="muted">No design conversations yet. Click "New ImmunoScope Session" to start.</p></div>';

  return `
    <div class="panel">
      <div class="form-header">
        <div>
          <h2>ImmunoScope</h2>
          <p class="muted">${tasks.length} conversation${tasks.length === 1 ? '' : 's'} — AI-assisted mutation design backed by MD evidence and literature</p>
        </div>
        <a class="button" href="#/copilot/new">New ImmunoScope Session</a>
      </div>
      <div class="grid cards" style="margin-top: 16px; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 16px;">
        ${taskCards}
      </div>
    </div>
  `;
}

function mountDesignCopilot() {
  document.querySelectorAll('[data-delete]').forEach(btn => {
    btn.addEventListener('click', async () => {
      const id = btn.getAttribute('data-delete');
      if (!confirm('Delete this design conversation? Drafts will be lost.')) return;
      try {
        await API.deleteDesignJob(id);
        window.App.navigate(false);
      } catch (e) {
        alert(`Failed: ${e.message}`);
      }
    });
  });
}

// ---------- New conversation page ----------

async function renderNewDesign() {
  const jobs = await API.listJobs().catch(() => []);
  const completedJobs = jobs.filter(j => j.status === 'completed');

  if (!completedJobs.length) {
    return `
      <div class="panel">
        <h2>New ImmunoScope Session</h2>
        <p class="muted">No completed analyses available.</p>
        <a class="button" href="#/new">Run an analysis first</a>
      </div>
    `;
  }

  const sourceOptions = completedJobs.map(j =>
    `<option value="${escapeHtml(j.id)}">${escapeHtml(j.name)} (${escapeHtml(j.id)})</option>`
  ).join('');

  const goalOptions = DESIGN_GOALS.map(g => `
    <label class="module-option">
      <input type="checkbox" name="design_goal" value="${escapeHtml(g.id)}">
      <span class="module-option-body">
        <span class="module-option-title">${escapeHtml(g.label)}</span>
        <span class="module-option-description">${escapeHtml(g.description)}</span>
      </span>
    </label>
  `).join('');

  // Task radio group. Default is TCR Affinity (matches the router default
  // and the historical web agent behavior).
  const taskOptions = DESIGN_TASKS.map((t, i) => `
    <label class="module-option">
      <input type="radio" name="task" value="${escapeHtml(t.id)}" ${i === 0 ? 'checked' : ''} required>
      <span class="module-option-body">
        <span class="module-option-title">${escapeHtml(t.label)}</span>
        <span class="module-option-description">${escapeHtml(t.description)}</span>
      </span>
    </label>
  `).join('');

  return `
    <form id="new-design-form" class="panel">
      <div class="form-header">
        <div>
          <h2>Start an ImmunoScope Session</h2>
          <p class="muted">Pick the source analysis and (optionally) set initial constraints. You'll be able to discuss design intent directly with the AI in the next step.</p>
        </div>
        <a href="#/copilot" style="color: var(--muted); text-decoration: none; font-size: 13px;">Cancel</a>
      </div>

      <label>
        Source analysis
        <select name="source_job_id" required>${sourceOptions}</select>
      </label>

      <label>
        Conversation name (optional)
        <input name="name" placeholder="e.g. cdr3_affinity_redesign">
      </label>

      <section class="module-picker">
        <div class="module-picker-head">
          <h3>Design task</h3>
        </div>
        <p class="muted" style="font-size: 12px; margin: 0 0 12px;">Pick which side of the interface you're mutating. This drives candidate filtering and the AI's framing — anchor-pocket chemistry is only surfaced for presentation tasks.</p>
        <div class="module-grid">${taskOptions}</div>
      </section>

      <section class="module-picker">
        <div class="module-picker-head">
          <h3>Design goals (optional hints)</h3>
        </div>
        <p class="muted" style="font-size: 12px; margin: 0 0 12px;">These nudge the AI's initial framing — you can change direction anytime in the conversation.</p>
        <div class="module-grid">${goalOptions}</div>
      </section>

      <section class="input-block">
        <h3>Initial constraints (optional)</h3>
        <div class="grid two-col">
          <label>
            Focus region
            <input name="focus_region" placeholder="e.g. CDR3_alpha, CDR3_beta">
          </label>
          <label>
            Must NOT mutate (comma-separated)
            <input name="must_exclude" placeholder="e.g. ASN96, CYS22">
          </label>
        </div>
        <label>
          Free-text notes for the AI
          <textarea name="notes" placeholder="e.g. Avoid introducing cysteines or prolines. Conservative mutations preferred." style="min-height: 60px;"></textarea>
        </label>
      </section>

      <button type="submit">Start conversation</button>
      <p id="new-design-status" class="muted"></p>
    </form>
  `;
}

function mountNewDesign() {
  const form = document.getElementById('new-design-form');
  const status = document.getElementById('new-design-status');

  form.addEventListener('submit', async event => {
    event.preventDefault();
    status.textContent = 'Creating conversation...';

    const formData = new FormData(form);
    const goals = formData.getAll('design_goal');
    const focusRegion = (formData.get('focus_region') || '').trim();
    const mustExclude = (formData.get('must_exclude') || '').split(',').map(s => s.trim()).filter(Boolean);
    const notes = (formData.get('notes') || '').trim();

    const constraints = {};
    if (focusRegion) constraints.focus_region = focusRegion;
    if (mustExclude.length) constraints.must_exclude = mustExclude;
    if (notes) constraints.notes = notes;

    const payload = {
      source_job_id: formData.get('source_job_id'),
      task: formData.get('task') || 'tcr_affinity',
      design_goals: goals,
      constraints,
      name: formData.get('name') || '',
      notes,
    };

    try {
      const result = await API.createDesignJob(payload);
      status.textContent = `Created ${result.id}. Opening conversation...`;
      setTimeout(() => { location.hash = `#/copilot/${result.id}`; }, 400);
    } catch (error) {
      status.textContent = `Failed: ${error.message}`;
    }
  });
}

// ---------- Three-column conversation workspace ----------

async function renderDesignDetail(taskId) {
  let task, drafts, context;
  try {
    const data = await API.getDesignJob(taskId);
    task = data.task;
    drafts = data.drafts || [];
    context = await API.getDesignContext(taskId);
  } catch (e) {
    return `<div class="panel"><h2>Error loading conversation</h2><p>${escapeHtml(e.message)}</p><a class="button" href="#/copilot">Back to list</a></div>`;
  }

  return `
    <div style="height: calc(100vh - 130px); display: grid; grid-template-columns: 280px minmax(0, 1fr) 320px; gap: 14px;">
      <!-- LEFT: source context -->
      <aside class="panel" style="overflow-y: auto; padding: 16px;">
        ${renderDesignContextPanel(task, context)}
      </aside>

      <!-- CENTER: conversation -->
      <section class="panel" style="display: flex; flex-direction: column; padding: 0; overflow: hidden;">
        <div style="padding: 14px 18px; border-bottom: 1px solid var(--line); display: flex; justify-content: space-between; align-items: center; gap: 12px;">
          <div style="flex: 1; min-width: 0;">
            <strong style="font-size: 14px;">${escapeHtml(task.name)}</strong>
            <p class="muted" style="margin: 2px 0 0; font-size: 11px;">
              <span style="display: inline-block; padding: 2px 6px; margin-right: 6px; border-radius: 4px; background: var(--accent-soft); color: var(--accent-strong); font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em;">${escapeHtml(DESIGN_TASK_LABEL[task.task] || 'tcr affinity')}</span>
              ImmunoScope · evidence-grounded mutation design
            </p>
          </div>
          <a href="#/copilot" style="color: var(--muted); text-decoration: none; font-size: 12px;">← All conversations</a>
        </div>
        <div id="design-chat-container" class="chat-container" style="flex: 1; overflow-y: auto;">
          ${task.task === 'peptide_presentation' ? `
          <div class="welcome-message" style="padding: 28px 20px;">
            <h3 style="margin: 0 0 12px;">Hi — I'm ImmunoScope (presentation track)</h3>
            <p style="font-size: 13px; line-height: 1.6;">I've loaded the MD evidence from <strong style="font-family: ui-monospace, monospace;">${escapeHtml(task.source_job_name)}</strong>. The anchor pocket chemistry for this allele is in my context — tell me what you want from the peptide (higher T½, preserve TCR recognition, shift selectivity) and I'll work through it with you.</p>
            <div class="example-questions" style="margin-top: 18px;">
              <button class="example-btn" data-example="Walk me through the peptide-HLA hotspots and explain which positions are anchors versus non-anchors.">Map peptide-HLA hotspots</button>
              <button class="example-btn" data-example="Suggest a conservative mutation at a non-anchor position that should improve presentation stability without breaking TCR recognition.">Improve T½ without breaking TCR</button>
              <button class="example-btn" data-example="If I had to mutate P2, which residues are tolerated on this allele and what is the expected impact?">Anchor-position alternatives</button>
            </div>
          </div>
          ` : `
          <div class="welcome-message" style="padding: 28px 20px;">
            <h3 style="margin: 0 0 12px;">Hi — I'm ImmunoScope</h3>
            <p style="font-size: 13px; line-height: 1.6;">I've loaded the MD evidence from <strong style="font-family: ui-monospace, monospace;">${escapeHtml(task.source_job_name)}</strong>. Tell me what you're trying to design — affinity, stability, specificity, or a more specific intent — and I'll work through it with you.</p>
            <div class="example-questions" style="margin-top: 18px;">
              <button class="example-btn" data-example="I want to enhance affinity by focusing on CDR3α. Walk me through the strongest hotspot candidates.">Enhance affinity via CDR3α</button>
              <button class="example-btn" data-example="What are the top 3 highest-confidence residues in this trajectory? Investigate them and save the strong picks.">Find top 3 high-confidence picks</button>
              <button class="example-btn" data-example="Show me hotspots, then pick the top candidate and explain why.">Show & analyze hotspots</button>
            </div>
          </div>
          `}
        </div>
        <div id="design-error-banner" class="error-banner" style="display: none; padding: 10px 18px;"></div>
        <div class="input-container" style="padding: 12px 16px;">
          <textarea id="design-chat-input" placeholder="Describe your design intent, ask for evidence, or refine recommendations..." rows="2"></textarea>
          <button id="design-send-btn">Send</button>
        </div>
      </section>

      <!-- RIGHT: drafts panel -->
      <aside class="panel" style="overflow: hidden; padding: 0; display: flex; flex-direction: column;">
        <div style="padding: 14px 16px; border-bottom: 1px solid var(--line); background: linear-gradient(180deg, rgba(236, 254, 255, 0.4), transparent);">
          <span style="font-size: 10px; font-weight: 800; color: var(--accent-strong); text-transform: uppercase; letter-spacing: 0.1em;">Draft recommendations</span>
          <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 4px;">
            <strong id="design-draft-count" style="font-size: 22px;">${drafts.length}</strong>
            <span class="muted" style="font-size: 11px;">saved by AI</span>
          </div>
          <a id="design-showcase-link" href="#/copilot/${escapeHtml(taskId)}/showcase" style="display: ${drafts.length ? 'inline-flex' : 'none'}; align-items: center; gap: 6px; margin-top: 10px; padding: 6px 12px; background: var(--accent); color: white; border-radius: 8px; text-decoration: none; font-size: 12px; font-weight: 700;">
            View Showcase →
          </a>
        </div>
        <div id="design-drafts-list" style="flex: 1; overflow-y: auto; padding: 12px;">
          ${renderDraftsList(drafts)}
        </div>
      </aside>
    </div>
  `;
}

function renderDesignContextPanel(task, context) {
  const hotspots = (context.top_hotspots || []).slice(0, 5);
  const hotspotsHtml = hotspots.length ? hotspots.map((h, i) => `
    <div style="display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px dashed var(--line); font-size: 12px;">
      <span style="font-family: ui-monospace, monospace; font-weight: 700;">#${i+1} ${escapeHtml(h.residue)}</span>
      <span class="muted" style="font-size: 11px;">${escapeHtml(h.region)} · RRCS ${h.rrcs}</span>
    </div>
  `).join('') : '<p class="muted" style="font-size: 12px;">No hotspots available</p>';

  const goalsHtml = (task.design_goals || []).length ? task.design_goals.map(g =>
    `<span style="display: inline-block; padding: 2px 8px; border-radius: 6px; background: var(--accent-soft); color: var(--accent-strong); font-size: 11px; font-weight: 700; margin-right: 4px; margin-bottom: 4px;">${escapeHtml(String(g).replace(/_/g, ' '))}</span>`
  ).join('') : '<span class="muted" style="font-size: 12px;">None set</span>';

  const constraints = task.constraints || {};
  const constraintsHtml = (constraints.must_exclude || constraints.focus_region || constraints.notes) ? `
    <div style="margin-top: 14px;">
      <span style="font-size: 10px; font-weight: 800; color: var(--muted); text-transform: uppercase; letter-spacing: 0.06em;">Constraints</span>
      <div style="margin-top: 6px; font-size: 12px; line-height: 1.5;">
        ${constraints.focus_region ? `<div><span class="muted">Focus:</span> <span style="font-family: ui-monospace, monospace;">${escapeHtml(constraints.focus_region)}</span></div>` : ''}
        ${constraints.must_exclude && constraints.must_exclude.length ? `<div><span class="muted">Exclude:</span> <span style="font-family: ui-monospace, monospace;">${escapeHtml(constraints.must_exclude.join(', '))}</span></div>` : ''}
        ${constraints.notes ? `<div style="margin-top: 4px; color: var(--text); font-style: italic;">${escapeHtml(constraints.notes)}</div>` : ''}
      </div>
    </div>
  ` : '';

  return `
    <div>
      <!-- Multi-agent design pipeline card (hydrated by initMultiAgentPanel) -->
      <div id="ma-pipeline-card" style="margin-bottom: 16px;"></div>

      <span style="font-size: 10px; font-weight: 800; color: var(--muted); text-transform: uppercase; letter-spacing: 0.1em;">Source system</span>
      <h3 style="margin: 4px 0 12px; font-family: ui-monospace, monospace; font-size: 14px; word-break: break-all;">${escapeHtml(task.source_job_name)}</h3>

      <div style="font-size: 12px; line-height: 1.6; padding-bottom: 12px; border-bottom: 1px solid var(--line);">
        ${context.peptide ? `<div><span class="muted">Peptide:</span> <span style="font-family: ui-monospace, monospace;">${escapeHtml(context.peptide)}</span></div>` : ''}
        ${context.hla ? `<div><span class="muted">HLA:</span> ${escapeHtml(context.hla)}</div>` : ''}
        ${context.tcr_genes ? `<div><span class="muted">TCR:</span> ${escapeHtml(context.tcr_genes)}</div>` : ''}
        <div><span class="muted">Frames:</span> ${context.n_frames || '-'}</div>
        <div><span class="muted">Duration:</span> ${context.duration_ns || '-'} ns</div>
        <div><span class="muted">Convergence:</span> <strong style="color: var(--accent-strong);">${escapeHtml(context.convergence || '-')}</strong></div>
        <div><span class="muted">RMSD:</span> ${context.rmsd || '-'} nm</div>
      </div>

      <div style="margin-top: 14px;">
        <span style="font-size: 10px; font-weight: 800; color: var(--muted); text-transform: uppercase; letter-spacing: 0.06em;">Design goals</span>
        <div style="margin-top: 6px;">${goalsHtml}</div>
      </div>

      ${constraintsHtml}

      <div style="margin-top: 18px;">
        <span style="font-size: 10px; font-weight: 800; color: var(--muted); text-transform: uppercase; letter-spacing: 0.06em;">Top hotspots</span>
        <div style="margin-top: 8px;">${hotspotsHtml}</div>
      </div>
    </div>
  `;
}

function renderDraftsList(drafts) {
  if (!drafts.length) {
    return `
      <div style="text-align: center; padding: 30px 12px; color: var(--muted);">
        <p style="font-size: 13px;">No recommendations saved yet.</p>
        <p style="font-size: 11px; margin-top: 8px;">When the AI commits to a recommendation, it'll appear here.</p>
      </div>
    `;
  }

  return drafts.map((d, idx) => {
    const pc = PRIORITY_COLORS[d.priority] || PRIORITY_COLORS.low;
    return `
      <div style="background: white; border: 1px solid var(--line); border-radius: 12px; padding: 12px; margin-bottom: 10px;">
        <div style="display: flex; justify-content: space-between; align-items: flex-start; gap: 8px; margin-bottom: 8px;">
          <div style="flex: 1;">
            <div style="font-family: ui-monospace, monospace; font-size: 14px; font-weight: 700;">${escapeHtml(chainRegionLabel(d) || d.residue || '?')}</div>
            <div class="muted" style="font-size: 10px; margin-top: 2px;">${escapeHtml(d.residue || '?')}</div>
          </div>
          <span style="padding: 2px 7px; border-radius: 999px; background: ${pc.bg}; color: ${pc.color}; font-size: 9px; font-weight: 800; text-transform: uppercase;">${escapeHtml(d.priority || 'low')}</span>
        </div>
        <div style="display: flex; align-items: center; gap: 6px; padding: 6px 8px; background: #f8fafc; border-radius: 8px; margin-bottom: 8px;">
          ${isPreserveDraft(d)
            ? `<span style="font-family: ui-monospace, monospace; font-size: 12px; color: var(--accent-strong); font-weight: 700;">Preserve ${escapeHtml(d.current_aa || '')}</span><span class="muted" style="font-size: 10px;">· do not mutate</span>`
            : `<span style="font-family: ui-monospace, monospace; color: var(--muted); font-size: 12px;">${escapeHtml(d.current_aa || '?')} →</span>
          ${(d.suggested_mutations || []).map(m => `<span style="font-family: ui-monospace, monospace; font-weight: 700; color: var(--accent-strong);">${escapeHtml(m)}</span>`).join(', ')}`}
        </div>
        ${d.rationale ? `<p style="font-size: 11px; line-height: 1.5; margin: 0 0 8px; color: var(--text);">${escapeHtml(d.rationale.length > 200 ? d.rationale.slice(0, 200) + '…' : d.rationale)}</p>` : ''}
        <div style="display: flex; justify-content: space-between; align-items: center; gap: 6px; padding-top: 6px; border-top: 1px solid var(--line);">
          <span class="muted" style="font-size: 10px;">conf: <strong style="color: var(--text); text-transform: uppercase;">${escapeHtml(d.confidence || '-')}</strong></span>
          <button data-draft-delete="${idx}" style="background: none; border: none; color: var(--muted); cursor: pointer; font-size: 11px; padding: 2px 6px;" title="Remove">×</button>
        </div>
      </div>
    `;
  }).join('');
}

// ---------- Multi-agent design pipeline (left-panel card + trace modal) ----------
// Decision A: the multi-agent pipeline is the recommender. This card drives a
// background run, streams stage progress (polling the status endpoint), drops
// the picks into the right-side drafts panel, and opens a full-trace modal.

let _maPollTimer = null;

function stopMaPolling() {
  if (_maPollTimer) { clearInterval(_maPollTimer); _maPollTimer = null; }
}

const MA_STAGE_LABEL = {
  start: 'Starting…',
  round1: 'Round 1 — parallel readers',
  round1_done: 'Round 1 complete',
  round2: 'Round 2 — recommender',
  critic: 'Round 2.5 — adversarial critic',
  revision: 'Revising recommendation',
  done: 'Done',
};

function _maEnsureSpinStyle() {
  if (document.getElementById('ma-spin-style')) return;
  const s = document.createElement('style');
  s.id = 'ma-spin-style';
  s.textContent = '@keyframes maspin{to{transform:rotate(360deg)}}'
    + ' .agent-process details>summary::-webkit-details-marker{display:none}'
    + ' .agent-process details>summary{list-style:none}'
    + ' .agent-process details[open]>summary>span:last-child::after{content:" ▾"}';
  document.head.appendChild(s);
}

async function initMultiAgentPanel(taskId) {
  stopMaPolling();
  _maEnsureSpinStyle();
  const card = document.getElementById('ma-pipeline-card');
  if (!card) return;
  let status = { status: 'idle' };
  try { status = await API.getDesignPipelineStatus(taskId); } catch (e) {}
  await renderMaCard(taskId, status);
  if (status.status === 'running') startMaPolling(taskId);
}

function startMaPolling(taskId) {
  stopMaPolling();
  _maPollTimer = setInterval(async () => {
    let status;
    try { status = await API.getDesignPipelineStatus(taskId); }
    catch (e) { return; }
    await renderMaCard(taskId, status);
    if (status.status !== 'running') {
      stopMaPolling();
      // Run finished — the picks were appended to drafts; refresh the panel.
      // (The chat picks up explain-mode context on its next WS connection.)
      if (window._designRefreshDrafts) { try { window._designRefreshDrafts(); } catch (e) {} }
      // Surface the recommended sites as rich cards in the conversation.
      if (window._designRenderPipelineCards) {
        try { window._designRenderPipelineCards(status.n_drafts_added || 0); } catch (e) {}
      }
    }
  }, 2000);
}

function maStagesStrip(status) {
  const events = (status && status.events) || [];
  const readersDone = events.filter(e => e.stage === 'round1' && e.status === 'done').length;
  const round2Done = events.some(e => e.stage === 'round2' && e.status === 'done');
  const criticSeen = events.some(e => e.stage === 'critic');
  const done = status.status && status.status !== 'running';
  const pill = (label, active, ok) => `
    <span style="display:inline-flex;align-items:center;gap:4px;padding:2px 8px;border-radius:999px;
      background:${ok ? '#dbe9e1' : active ? 'var(--accent-soft)' : '#f1f0ec'};
      color:${ok ? '#255b4b' : active ? 'var(--accent-strong)' : '#9a958c'};
      font-size:10px;font-weight:800;letter-spacing:0.03em;">${ok ? '✓' : '•'} ${label}</span>`;
  return `
    <div style="display:flex;flex-wrap:wrap;gap:4px;margin-top:8px;">
      ${pill(`readers ${readersDone}/4`, readersDone > 0 && readersDone < 4, readersDone >= 4)}
      ${pill('recommend', !round2Done && readersDone >= 4, round2Done)}
      ${pill('critic', criticSeen && !done, criticSeen && done)}
      ${pill('done', false, done)}
    </div>`;
}

function maRunningHtml(status) {
  const progress = Math.max(3, Math.min(100, status.progress || 5));
  const stageLabel = MA_STAGE_LABEL[status.stage] || status.stage || 'running';
  return `
    <div style="border:1px solid var(--accent); border-radius:12px; padding:12px; background:linear-gradient(180deg, rgba(236,254,255,0.5), transparent);">
      <div style="display:flex;align-items:center;gap:8px;">
        <span style="width:12px;height:12px;border:2px solid var(--accent);border-top-color:transparent;border-radius:50%;display:inline-block;animation:maspin 0.8s linear infinite;"></span>
        <strong style="font-size:12px;">Multi-agent design running</strong>
      </div>
      <div style="margin-top:8px;height:6px;background:#e8f3f1;border-radius:999px;overflow:hidden;">
        <div style="height:100%;width:${progress}%;background:var(--accent);transition:width 0.4s;"></div>
      </div>
      <div class="muted" style="font-size:11px;margin-top:6px;">${escapeHtml(stageLabel)} · ${progress}%</div>
      ${maStagesStrip(status)}
    </div>`;
}

function maSummaryHtml(status, art) {
  if (!art) {
    const err = status.error ? `<div style="font-size:11px;color:#9b3a2a;margin-top:6px;">${escapeHtml(status.error)}</div>` : '';
    return `<div style="font-size:12px;margin-top:8px;">Status: <strong>${escapeHtml(status.status || 'unknown')}</strong></div>${err}`;
  }
  const rec = (art.recommendation && art.recommendation.structured_payload) || {};
  const recs = rec.recommendations || [];
  const verdict = (art.critique && art.critique.verdict) || '—';
  const det = (art.deterministic_findings || []).length;
  const vColor = verdict === 'approve' ? '#255b4b' : verdict === 'revise' ? '#a36b00' : '#6f746d';
  return `
    <div style="margin-top:8px;display:grid;grid-template-columns:1fr 1fr;gap:6px;font-size:12px;">
      <div style="padding:6px 8px;background:white;border:1px solid var(--line);border-radius:8px;">
        <div class="muted" style="font-size:9px;text-transform:uppercase;letter-spacing:0.06em;">Site picks</div>
        <strong style="font-size:15px;">${recs.length}</strong>
      </div>
      <div style="padding:6px 8px;background:white;border:1px solid var(--line);border-radius:8px;">
        <div class="muted" style="font-size:9px;text-transform:uppercase;letter-spacing:0.06em;">Critic</div>
        <strong style="font-size:13px;color:${vColor};text-transform:capitalize;">${escapeHtml(verdict)}</strong>
      </div>
    </div>
    ${det ? `<div class="muted" style="font-size:10px;margin-top:6px;">${det} machine-confirmed finding(s)</div>` : ''}`;
}

async function renderMaCard(taskId, status) {
  const card = document.getElementById('ma-pipeline-card');
  if (!card) return;
  const st = (status && status.status) || 'idle';
  if (st === 'running') { card.innerHTML = maRunningHtml(status); return; }

  let summaryHtml = '';
  let runLabel = 'Run multi-agent design';
  if (st === 'completed' || st === 'completed_with_errors' || st === 'failed') {
    runLabel = 'Re-run pipeline';
    let art = null;
    try { const r = await API.getDesignPipeline(taskId); if (r && r.available) art = r.artifact; } catch (e) {}
    summaryHtml = maSummaryHtml(status, art);
  }
  const tip = st === 'idle'
    ? `<p class="muted" style="font-size:11px;line-height:1.5;margin:8px 0 0;">Four parallel readers → recommender → adversarial critic. Produces the authoritative site picks (right panel) plus a full trace. Takes a few minutes.</p>`
    : '';
  card.innerHTML = `
    <div style="border:1px solid var(--line);border-radius:12px;padding:12px;background:#fbfdfc;">
      <div style="display:flex;align-items:center;justify-content:space-between;gap:8px;">
        <span style="font-size:10px;font-weight:800;color:var(--accent-strong);text-transform:uppercase;letter-spacing:0.08em;">⚙ Multi-agent design</span>
        ${st !== 'idle' ? `<button data-ma-trace style="background:none;border:none;color:var(--accent-strong);font-size:11px;font-weight:700;cursor:pointer;padding:0;">View full trace →</button>` : ''}
      </div>
      ${summaryHtml}
      ${tip}
      <button data-ma-run style="margin-top:10px;width:100%;padding:8px 12px;border:none;border-radius:8px;background:var(--accent);color:white;font-size:12px;font-weight:700;cursor:pointer;">${runLabel}</button>
    </div>`;
  bindMaButtons(taskId);
}

function bindMaButtons(taskId) {
  const card = document.getElementById('ma-pipeline-card');
  if (!card) return;
  const runBtn = card.querySelector('[data-ma-run]');
  if (runBtn) {
    runBtn.addEventListener('click', async () => {
      runBtn.disabled = true; runBtn.textContent = 'Starting…';
      try {
        await API.runDesignPipeline(taskId);
        await renderMaCard(taskId, { status: 'running', progress: 5, stage: 'start', events: [] });
        startMaPolling(taskId);
      } catch (e) {
        runBtn.disabled = false; runBtn.textContent = 'Run multi-agent design';
        alert(`Failed to start pipeline: ${e.message}`);
      }
    });
  }
  const traceBtn = card.querySelector('[data-ma-trace]');
  if (traceBtn) traceBtn.addEventListener('click', () => openMaTraceModal(taskId));
}

function _maPre(txt) {
  return `<div style="white-space:pre-wrap;font-size:12px;line-height:1.55;background:#fbfdfc;border:1px solid var(--line);border-radius:8px;padding:10px;">${escapeHtml(txt || '')}</div>`;
}

function maTraceBody(art) {
  const sec = (title, inner) => `<div style="margin-bottom:18px;"><div style="font-size:11px;font-weight:800;text-transform:uppercase;letter-spacing:0.08em;color:var(--muted);margin-bottom:8px;">${title}</div>${inner}</div>`;

  const statusLine = `<div style="display:flex;gap:8px;flex-wrap:wrap;margin-bottom:16px;">
    <span style="padding:3px 10px;border-radius:999px;background:${art.succeeded ? '#dbe9e1' : '#f3e8cc'};color:${art.succeeded ? '#255b4b' : '#a36b00'};font-size:11px;font-weight:800;">${art.succeeded ? 'succeeded' : 'completed with errors'}</span>
    ${art.critique && art.critique.verdict ? `<span style="padding:3px 10px;border-radius:999px;background:#eaf0ff;color:#3a4a7a;font-size:11px;font-weight:800;">critic: ${escapeHtml(art.critique.verdict)}</span>` : ''}
  </div>`;

  const rec = art.recommendation || {};
  const recPayload = rec.structured_payload || {};
  const recs = recPayload.recommendations || [];
  let recTable = '';
  if (recs.length) {
    recTable = `<table style="width:100%;border-collapse:collapse;font-size:12px;margin-bottom:10px;">
      <thead><tr style="text-align:left;color:var(--muted);">
        <th style="padding:6px;border-bottom:1px solid var(--line);">Residue</th>
        <th style="padding:6px;border-bottom:1px solid var(--line);">Role</th>
        <th style="padding:6px;border-bottom:1px solid var(--line);">Priority</th>
        <th style="padding:6px;border-bottom:1px solid var(--line);">Direction</th>
        <th style="padding:6px;border-bottom:1px solid var(--line);">Candidate AAs</th>
        <th style="padding:6px;border-bottom:1px solid var(--line);">Rationale</th>
      </tr></thead><tbody>
      ${recs.map(r => { const h = (r.handoff_to_physchem || {}); return `<tr>
        <td style="padding:6px;border-bottom:1px solid #f0ebe0;font-family:ui-monospace,monospace;font-weight:700;">${escapeHtml(r.residue || '?')}</td>
        <td style="padding:6px;border-bottom:1px solid #f0ebe0;">${escapeHtml(r.role || '-')}</td>
        <td style="padding:6px;border-bottom:1px solid #f0ebe0;">${escapeHtml(r.priority || '-')}</td>
        <td style="padding:6px;border-bottom:1px solid #f0ebe0;font-size:11px;">${escapeHtml(h.design_hint_class || '-')}</td>
        <td style="padding:6px;border-bottom:1px solid #f0ebe0;font-family:ui-monospace,monospace;font-size:11px;">${escapeHtml((h.candidate_residues || []).join(', ') || '-')}</td>
        <td style="padding:6px;border-bottom:1px solid #f0ebe0;font-size:11px;">${escapeHtml(r.rationale_brief || '')}</td>
      </tr>`; }).join('')}
      </tbody></table>`;
  }
  const skipped = recPayload.skipped_candidates || [];
  const skipHtml = skipped.length ? `<div class="muted" style="font-size:11px;margin-top:6px;">Skipped: ${skipped.map(s => `${escapeHtml(s.residue || '?')} (${escapeHtml(s.reason || '')})`).join('; ')}</div>` : '';
  const recSection = sec('Recommendation (authoritative)', recTable + (rec.narrative ? _maPre(rec.narrative) : '') + skipHtml);

  let criticInner;
  if (art.critique) {
    const findings = art.critique.findings || [];
    criticInner = (art.critique.narrative ? _maPre(art.critique.narrative) : '') +
      (findings.length ? `<ul style="margin:8px 0 0;padding-left:18px;font-size:12px;line-height:1.5;">${findings.map(f => `<li>[${escapeHtml(String(f.severity || '?'))}] ${escapeHtml(f.issue || f.detail || f.description || JSON.stringify(f))}</li>`).join('')}</ul>` : '<p class="muted" style="font-size:12px;">No findings.</p>');
  } else { criticInner = '<p class="muted" style="font-size:12px;">Critic did not run.</p>'; }
  const criticSection = sec('Design critic (round 2.5)', criticInner);

  const det = art.deterministic_findings || [];
  const detSection = det.length ? sec('Machine-confirmed findings (model-free)',
    `<ul style="margin:0;padding-left:18px;font-size:12px;line-height:1.5;">${det.map(d => `<li>[${escapeHtml(d.severity || '?')}] ${escapeHtml(d.issue_type || '')} @ ${escapeHtml(d.target || '')}: ${escapeHtml(d.detail || '')}</li>`).join('')}</ul>`) : '';

  const audit = art.numeric_audit_summary ? sec('Numeric-fidelity audit (round-1)', _maPre(art.numeric_audit_summary)) : '';
  const cand = art.candidate_table ? sec('Deterministic candidate table', _maPre(art.candidate_table)) : '';

  const briefs = art.briefs || {};
  const briefOrder = [['bio', 'Biology'], ['interaction', 'Interaction'], ['conformation', 'Conformation'], ['interface_exposure', 'Interface & exposure']];
  const briefsInner = briefOrder.map(([k, label]) => {
    const b = briefs[k];
    if (!b) return `<details style="margin-bottom:6px;"><summary style="cursor:pointer;font-weight:700;font-size:12px;color:#9b3a2a;">${label} — failed</summary></details>`;
    return `<details style="margin-bottom:6px;"><summary style="cursor:pointer;font-weight:700;font-size:12px;">${label}</summary>
      <div style="margin-top:6px;">${_maPre(b.narrative)}</div>
      ${b.numeric_audit ? `<div class="muted" style="font-size:10px;margin-top:4px;">audit: ${escapeHtml(b.numeric_audit)}</div>` : ''}
    </details>`;
  }).join('');
  const briefsSection = sec('Round-1 reader briefs', briefsInner);

  return statusLine + recSection + criticSection + detSection + audit + cand + briefsSection;
}

async function openMaTraceModal(taskId) {
  let art = null;
  try { const r = await API.getDesignPipeline(taskId); if (r && r.available) art = r.artifact; } catch (e) {}
  const overlay = document.createElement('div');
  overlay.style.cssText = 'position:fixed;inset:0;background:rgba(15,23,20,0.55);z-index:9999;display:flex;align-items:center;justify-content:center;padding:24px;';
  overlay.addEventListener('click', (e) => { if (e.target === overlay) overlay.remove(); });
  overlay.innerHTML = `
    <div style="background:#fff;border-radius:16px;max-width:980px;width:100%;max-height:90vh;overflow-y:auto;box-shadow:0 20px 60px rgba(0,0,0,0.3);">
      <div style="position:sticky;top:0;background:#fff;padding:16px 22px;border-bottom:1px solid var(--line);display:flex;justify-content:space-between;align-items:center;z-index:1;">
        <strong style="font-size:15px;">Multi-agent design trace</strong>
        <button data-ma-close style="background:none;border:none;font-size:22px;cursor:pointer;color:var(--muted);line-height:1;">×</button>
      </div>
      <div style="padding:18px 22px;">${art ? maTraceBody(art) : '<p class="muted">No trace artifact available yet — run the pipeline first.</p>'}</div>
    </div>`;
  document.body.appendChild(overlay);
  overlay.querySelector('[data-ma-close]').addEventListener('click', () => overlay.remove());
}

// Mount design conversation (called after renderDesignDetail)
function mountDesignDetail(taskId) {
  const chatContainer = document.getElementById('design-chat-container');
  const inputEl = document.getElementById('design-chat-input');
  const sendBtn = document.getElementById('design-send-btn');
  const draftsList = document.getElementById('design-drafts-list');
  const draftCount = document.getElementById('design-draft-count');
  const errorBanner = document.getElementById('design-error-banner');

  if (!chatContainer) return;

  // Open WebSocket with design_task_id parameter so the backend injects
  // design context into the session's system prompt.
  const wsProto = location.protocol === 'https:' ? 'wss:' : 'ws:';
  const sessionId = `design_${taskId}_${Date.now()}`;
  const wsUrl = `${wsProto}//${location.host}/api/agent/ws?session_id=${encodeURIComponent(sessionId)}&design_task_id=${encodeURIComponent(taskId)}`;
  let ws = null;
  let currentAgentMsg = null;
  let agentBuffer = '';
  // Option A: per-turn collapsible "Working" strip. Tool calls + intermediate
  // narration fold into one muted, collapsed block so the reply reads clean.
  let currentProcessBlock = null;
  let processStepCount = 0;

  function appendMessage(role, content) {
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    const wrap = document.createElement('div');
    wrap.className = `message ${role}`;
    const avatar = role === 'user' ? 'U' : 'AI';
    wrap.innerHTML = `
      <div class="message-avatar">${avatar}</div>
      <div class="message-content"></div>
    `;
    wrap.querySelector('.message-content').textContent = content;
    chatContainer.appendChild(wrap);
    chatContainer.scrollTop = chatContainer.scrollHeight;
    return wrap.querySelector('.message-content');
  }

  function appendToolCall(toolName, status) {
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    const wrap = document.createElement('div');
    wrap.style.cssText = 'padding: 4px 16px; display: flex; gap: 12px;';
    wrap.innerHTML = `
      <div class="message-avatar" style="background: var(--accent-soft); color: var(--accent); width: 36px; height: 36px;">⚙</div>
      <div class="tool-call ${status}" style="font-family: ui-monospace, monospace;">
        <span>${escapeHtml(toolName)}</span>
      </div>
    `;
    chatContainer.appendChild(wrap);
    chatContainer.scrollTop = chatContainer.scrollHeight;
    return wrap;
  }

  // ---- Agent "process" fold (option A) ----------------------------------
  // Collect a turn's tool calls + intermediate narration into ONE collapsible
  // muted strip, so the user-facing reply (the final text bubble) reads clean.
  function getProcessBlock() {
    if (currentProcessBlock && chatContainer.contains(currentProcessBlock)) return currentProcessBlock;
    if (typeof _maEnsureSpinStyle === 'function') _maEnsureSpinStyle();
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    const wrap = document.createElement('div');
    wrap.className = 'agent-process';
    wrap.style.cssText = 'padding: 2px 16px 2px 60px;';
    wrap.innerHTML = `
      <details style="border:1px solid var(--line); border-radius:8px; background:#fafafa;">
        <summary style="cursor:pointer; list-style:none; display:flex; align-items:center; gap:7px; padding:6px 10px; font-size:12px; color:var(--muted); user-select:none;">
          <span class="ap-spin" style="display:inline-block; width:12px; height:12px; border:2px solid var(--muted); border-top-color:transparent; border-radius:50%; animation:maspin 0.8s linear infinite;"></span>
          <span>Working · <span class="ap-count">0</span> steps</span>
          <span style="margin-left:auto; font-size:10px; opacity:0.7;">expand</span>
        </summary>
        <div class="ap-body" style="padding:6px 12px 8px 28px; font-family: ui-monospace, monospace; font-size:11px; color:#6f746d; line-height:1.7; border-top:1px solid var(--line); white-space:pre-wrap; word-break:break-word;"></div>
      </details>`;
    chatContainer.appendChild(wrap);
    currentProcessBlock = wrap;
    processStepCount = 0;
    return wrap;
  }

  function addProcessStep(label) {
    const block = getProcessBlock();
    const line = document.createElement('div');
    line.textContent = label;
    block.querySelector('.ap-body').appendChild(line);
    processStepCount += 1;
    const cnt = block.querySelector('.ap-count');
    if (cnt) cnt.textContent = String(processStepCount);
    chatContainer.scrollTop = chatContainer.scrollHeight;
    return line;
  }

  // Move an in-progress narration bubble ("let me check…") into the process
  // fold, so only the final answer survives as a clean reply bubble.
  function foldNarrationIntoProcess() {
    if (!currentAgentMsg) return;
    const txt = (agentBuffer || '').trim();
    const wrap = currentAgentMsg.closest('.message');
    if (txt) addProcessStep('💬 ' + txt);
    if (wrap) wrap.remove();
    currentAgentMsg = null;
    agentBuffer = '';
  }

  // End-of-turn: stop the spinner, mark done, and close the block so the next
  // turn opens a fresh one.
  function finalizeProcessBlock() {
    if (!currentProcessBlock) return;
    const spin = currentProcessBlock.querySelector('.ap-spin');
    if (spin) {
      spin.style.animation = 'none';
      spin.style.border = 'none';
      spin.style.width = 'auto';
      spin.style.height = 'auto';
      spin.style.color = '#255b4b';
      spin.style.fontWeight = '800';
      spin.textContent = '✓';
    }
    currentProcessBlock = null;
    processStepCount = 0;
  }

  function appendMutationProfileCard(rec) {
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    const wrap = document.createElement('div');
    wrap.className = 'mutation-profile-card-wrap';
    wrap.style.cssText = 'padding: 8px 16px;';
    wrap.innerHTML = renderMutationProfileCard(rec);
    chatContainer.appendChild(wrap);
    // Draw molecules via SmilesDrawer
    renderPendingMolecules(wrap);

    // Initialize inline 3D viewer + bind "Expand" button.
    // The PDB endpoint is keyed by design task id (not source job id).
    wrap.querySelectorAll('[data-inline-viewer]').forEach(viewerEl => {
      mountInline3DViewer(viewerEl, rec, taskId);
    });
    wrap.querySelectorAll('.btn-view-3d').forEach(btn => {
      btn.addEventListener('click', () => open3DViewer(rec, taskId));
    });

    // Populate backbone flexibility slot (φ/ψ from cached context)
    getDesignContextCached(taskId).then(ctx => {
      const slot = wrap.querySelector('[data-backbone-flex-slot]');
      if (!slot) return;
      const resIndex = (ctx && ctx.dihedrals && ctx.dihedrals.residue_index) || {};
      const info = resIndex[rec.residue];
      const dockingAngles = (ctx && ctx.docking_angles) || {};
      if (!info && !dockingAngles.crossing_mean) return;
      slot.style.display = 'block';
      const spread = info ? info.angular_spread_deg : null;
      const flexLabel = spread === null ? '-'
        : spread > 30 ? 'highly flexible'
        : spread > 15 ? 'moderately flexible'
        : 'rigid';
      const flexColor = spread === null ? '#6f746d'
        : spread > 30 ? '#255b4b'
        : spread > 15 ? '#a36b00'
        : '#9b3a2a';
      slot.innerHTML = `
        <div style="padding: 12px 18px; border-bottom: 1px solid #f0ebe0;">
          <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d; margin-bottom: 8px;">Backbone geometry context</div>
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px;">
            ${info ? `
              <div style="padding: 10px; background: #fffdfa; border: 1px solid #ddd5c7; border-radius: 10px;">
                <div style="font-size: 10px; color: #6f746d; font-weight: 700; letter-spacing: 0.06em; text-transform: uppercase;">Backbone (φ/ψ)</div>
                <div style="margin-top: 4px; font-family: ui-monospace, monospace; font-size: 13px;">${info.angular_spread_deg.toFixed(1)}° spread</div>
                <div style="margin-top: 2px; font-size: 11px; color: ${flexColor}; font-weight: 700;">${flexLabel}</div>
                <div style="margin-top: 4px; font-size: 10px; color: #6f746d;">${info.dominant_region} · φ̄=${info.phi_mean.toFixed(0)}° ψ̄=${info.psi_mean.toFixed(0)}°</div>
              </div>
            ` : '<div></div>'}
            ${dockingAngles.crossing_mean !== null && dockingAngles.crossing_mean !== undefined ? `
              <div style="padding: 10px; background: #fffdfa; border: 1px solid #ddd5c7; border-radius: 10px;">
                <div style="font-size: 10px; color: #6f746d; font-weight: 700; letter-spacing: 0.06em; text-transform: uppercase;">TCR-pMHC docking</div>
                <div style="margin-top: 4px; font-family: ui-monospace, monospace; font-size: 13px;">Crossing ${dockingAngles.crossing_mean.toFixed(1)}°</div>
                <div style="margin-top: 2px; font-size: 11px; color: ${dockingAngles.crossing_std < 5 ? '#255b4b' : '#a36b00'}; font-weight: 700;">σ=${(dockingAngles.crossing_std||0).toFixed(1)}° ${dockingAngles.crossing_std < 5 ? 'stable' : 'flexible'}</div>
                ${dockingAngles.incident_mean !== null && dockingAngles.incident_mean !== undefined ? `<div style="margin-top: 4px; font-size: 10px; color: #6f746d;">Incident ${dockingAngles.incident_mean.toFixed(1)}°</div>` : ''}
              </div>
            ` : '<div></div>'}
          </div>
        </div>
      `;
    }).catch(err => console.warn('Backbone flex hydration failed:', err));

    chatContainer.scrollTop = chatContainer.scrollHeight;
    return wrap;
  }

  // Lightweight residue viewer card — rendered when Agent calls
  // visualize_residue tool. Smaller and focused on the 3D view, not a full
  // mutation profile.
  function appendResidueViewerCard(payload) {
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    const wrap = document.createElement('div');
    wrap.className = 'residue-viewer-card-wrap';
    wrap.style.cssText = 'padding: 8px 16px;';

    // Build a "fake" rec-like object so we can reuse mountInline3DViewer +
    // extractContactsFromEvidence with minimal changes.
    const rec = {
      residue: payload.residue || '?',
      chain: payload.chain || '',
      region: payload.region || '',
      suggested_mutations: payload.suggested_mutations || [],
      // Convert contact_partners to supporting_evidence format that
      // extractContactsFromEvidence understands (md_data with "RESIDUE (chain)" patterns)
      supporting_evidence: (payload.contact_partners || []).map(p => ({
        type: 'md_data',
        data: `${p.residue || '?'} (${p.chain || 'peptide'})`,
      })),
      // Also keep interaction_partners as structured field
      interaction_partners: payload.contact_partners || [],
    };

    const mutationsHTML = (rec.suggested_mutations || []).length
      ? `<div style="font-size: 11px; color: #6f746d; margin-top: 4px;">Considered mutations: ${rec.suggested_mutations.map(m => `<span style="display: inline-block; padding: 1px 6px; border-radius: 4px; background: var(--accent-soft); color: var(--accent-strong); font-family: ui-monospace, monospace; font-weight: 700; margin: 0 2px;">${escapeHtml(m)}</span>`).join('')}</div>`
      : '';

    wrap.innerHTML = `
      <div style="background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%); border: 1px solid #ddd5c7; border-radius: 16px; box-shadow: 0 8px 20px rgba(29, 37, 32, 0.06); overflow: hidden; max-width: 720px; margin-left: 52px; font-family: 'Noto Sans', 'PingFang SC', sans-serif; color: #1e2a24;">
        <div style="padding: 10px 16px; border-bottom: 1px solid #f0ebe0; background: rgba(245, 241, 232, 0.45); display: flex; align-items: center; gap: 10px;">
          <span style="font-size: 9px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d;">3D inspection</span>
          <span style="font-family: ui-monospace, monospace; font-size: 14px; font-weight: 800; color: #1e2a24;">${escapeHtml(rec.residue)}</span>
          ${rec.region ? `<span style="font-size: 11px; color: #6f746d;">${escapeHtml(rec.chain || '')} · ${escapeHtml(rec.region || '')}</span>` : ''}
          <button class="btn-expand-3d" style="margin-left: auto; padding: 3px 10px; border-radius: 999px; background: #fffdfa; border: 1px solid #ddd5c7; color: #255b4b; font-size: 10px; font-weight: 800; cursor: pointer; font-family: 'Noto Sans', sans-serif; letter-spacing: 0.06em;">Expand →</button>
        </div>
        <div style="padding: 12px 16px;">
          ${renderInline3DViewer(rec)}
          <div style="display: flex; justify-content: center; gap: 18px; margin-top: 8px; font-size: 10px; color: #6f746d;">
            <span style="display: inline-flex; align-items: center; gap: 5px;"><span style="display: inline-block; width: 10px; height: 10px; border-radius: 50%; background: #dc2626;"></span>Target residue</span>
            <span style="display: inline-flex; align-items: center; gap: 5px;"><span style="display: inline-block; width: 10px; height: 10px; border-radius: 50%; background: #f59e0b;"></span>Contact partner</span>
            <span style="display: inline-flex; align-items: center; gap: 5px;"><span style="display: inline-block; width: 14px; border-top: 1.5px dashed #fde047;"></span>H-bond</span>
          </div>
          ${mutationsHTML}
          ${payload.note ? `<p style="margin: 10px 0 0; padding: 8px 12px; background: var(--accent-soft); color: var(--accent-strong); border-radius: 8px; font-size: 12px; line-height: 1.5; font-style: italic;">${escapeHtml(payload.note)}</p>` : ''}
        </div>
      </div>
    `;
    chatContainer.appendChild(wrap);

    // Initialize the 3D viewer + bind expand button.
    // PDB endpoint is keyed by design task id.
    wrap.querySelectorAll('[data-inline-viewer]').forEach(viewerEl => {
      mountInline3DViewer(viewerEl, rec, taskId);
    });
    wrap.querySelectorAll('.btn-expand-3d').forEach(btn => {
      btn.addEventListener('click', () => open3DViewer(rec, taskId));
    });

    chatContainer.scrollTop = chatContainer.scrollHeight;
    return wrap;
  }

  // Comparison-ready card — shown when compare_systems tool finishes.
  // Compact "you have a new comparison" card with quick links.
  function appendComparisonReadyCard(payload) {
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    const wrap = document.createElement('div');
    wrap.className = 'comparison-ready-wrap';
    wrap.style.cssText = 'padding: 8px 16px;';
    wrap.innerHTML = `
      <div style="background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%); border: 1px solid #ddd5c7; border-radius: 16px; box-shadow: 0 10px 24px rgba(29, 37, 32, 0.06); overflow: hidden; max-width: 720px; margin-left: 52px; font-family: 'Noto Sans', 'PingFang SC', sans-serif; color: #1e2a24;">
        <div style="padding: 14px 18px; background: linear-gradient(135deg, #1c473c 0%, #2d6f5c 100%); color: #f7f5ef;">
          <span style="font-size: 9px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; opacity: 0.85;">Comparison ready</span>
          <div style="margin-top: 4px; font-size: 16px; font-weight: 800;">${escapeHtml(payload.name || 'Comparison')}</div>
          <div style="margin-top: 4px; font-size: 12px; opacity: 0.85;">
            <strong>${escapeHtml(payload.label_a || 'A')}</strong> vs <strong>${escapeHtml(payload.label_b || 'B')}</strong>
          </div>
        </div>
        <div style="padding: 14px 18px; display: flex; gap: 8px; flex-wrap: wrap; align-items: center;">
          <a href="#/jobs/${escapeHtml(payload.compare_job_id)}" style="flex: 1; padding: 10px 14px; border-radius: 10px; background: #1c473c; color: #fffdf8; text-decoration: none; font-weight: 800; text-align: center; font-size: 13px;">Open Comparison →</a>
          <button class="btn-show-card" data-job-id="${escapeHtml(payload.compare_job_id)}" style="padding: 10px 14px; border-radius: 10px; background: #fffdfa; border: 1px solid #ddd5c7; color: #255b4b; font-weight: 800; cursor: pointer; font-size: 13px;">Show summary here</button>
        </div>
      </div>
    `;
    chatContainer.appendChild(wrap);
    // "Show summary here" → fetch + render comparison card inline
    wrap.querySelectorAll('.btn-show-card').forEach(btn => {
      btn.addEventListener('click', () => {
        const id = btn.getAttribute('data-job-id');
        appendComparisonCard({ compare_job_id: id, focus: 'overview' });
      });
    });
    chatContainer.scrollTop = chatContainer.scrollHeight;
    return wrap;
  }

  // Inline comparison summary card — rendered by show_comparison_card tool
  // or via the "Show summary here" button. Pulls real data from the API.
  function appendComparisonCard(payload) {
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    const wrap = document.createElement('div');
    wrap.className = 'comparison-card-wrap';
    wrap.style.cssText = 'padding: 8px 16px;';
    wrap.innerHTML = `
      <div style="background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%); border: 1px solid #ddd5c7; border-radius: 16px; box-shadow: 0 10px 24px rgba(29, 37, 32, 0.06); overflow: hidden; max-width: 720px; margin-left: 52px; font-family: 'Noto Sans', 'PingFang SC', sans-serif; color: #1e2a24; min-height: 100px;">
        <div style="padding: 24px; text-align: center; color: #6f746d; font-size: 13px;">
          Loading comparison summary...
        </div>
      </div>
    `;
    chatContainer.appendChild(wrap);
    chatContainer.scrollTop = chatContainer.scrollHeight;

    const compareJobId = payload.compare_job_id;
    Promise.all([
      API.getCompareSummary(compareJobId),
      API.listComparePlots(compareJobId),
    ]).then(([summaryResp, plotsResp]) => {
      if (summaryResp.status !== 'ready') {
        wrap.innerHTML = `
          <div style="background: #fffdfa; border: 1px solid #ddd5c7; border-radius: 16px; padding: 18px; max-width: 720px; margin-left: 52px;">
            <div style="font-size: 11px; font-weight: 800; color: #6f746d; letter-spacing: 0.1em; text-transform: uppercase;">Comparison not ready</div>
            <p style="margin: 8px 0 0; font-size: 13px; color: #1e2a24;">${escapeHtml(summaryResp.message || 'Artifacts not generated yet.')}</p>
          </div>
        `;
        return;
      }
      const summary = summaryResp.summary || {};
      const caseA = summary.case_a || {};
      const caseB = summary.case_b || {};
      const takeaways = summary.takeaways || [];
      const comparability = summary.comparability || {};
      const plots = plotsResp.plots || [];

      // Pick 4 key plots for the inline card (smaller, more compact)
      const keyPlots = plots.slice(0, 4);

      const compStatus = (comparability.status || '').toLowerCase();
      const compBadge = compStatus === 'comparable'
        ? { bg: '#dbe9e1', color: '#255b4b', label: 'Comparable' }
        : compStatus === 'caution'
          ? { bg: '#f3e8cc', color: '#a36b00', label: 'Caution' }
          : { bg: '#f0ddd5', color: '#9b3a2a', label: comparability.status || 'Check' };

      wrap.innerHTML = `
        <div style="background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%); border: 1px solid #ddd5c7; border-radius: 16px; box-shadow: 0 10px 24px rgba(29, 37, 32, 0.06); overflow: hidden; max-width: 720px; margin-left: 52px; font-family: 'Noto Sans', 'PingFang SC', sans-serif; color: #1e2a24;">
          <!-- Header -->
          <div style="padding: 12px 18px; background: linear-gradient(180deg, rgba(245, 241, 232, 0.55), transparent); border-bottom: 1px solid #ddd5c7;">
            <div style="display: flex; align-items: center; gap: 8px; flex-wrap: wrap;">
              <span style="font-size: 9px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d;">Comparison summary</span>
              <span style="padding: 2px 8px; border-radius: 999px; background: ${compBadge.bg}; color: ${compBadge.color}; font-size: 9px; font-weight: 800; letter-spacing: 0.06em;">${escapeHtml(compBadge.label)}</span>
            </div>
            <div style="font-size: 14px; font-weight: 700; margin-top: 6px;">
              <strong style="color: #255b4b;">${escapeHtml(caseA.label || 'A')}</strong>
              <span style="color: #6f746d; margin: 0 6px;">vs</span>
              <strong style="color: #a36b00;">${escapeHtml(caseB.label || 'B')}</strong>
            </div>
          </div>

          <!-- Takeaways -->
          ${takeaways.length ? `
            <div style="padding: 12px 18px; border-bottom: 1px solid #f0ebe0;">
              <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #6f746d; margin-bottom: 6px;">Key takeaways</div>
              <ul style="margin: 0; padding-left: 18px; font-size: 13px; line-height: 1.6;">
                ${takeaways.map(t => `<li>${escapeHtml(t)}</li>`).join('')}
              </ul>
            </div>
          ` : ''}

          <!-- Plot grid (4 small thumbnails) -->
          ${keyPlots.length ? `
            <div style="padding: 12px 18px; border-bottom: 1px solid #f0ebe0; display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px;">
              ${keyPlots.map(p => `
                <div style="text-align: center;">
                  <div style="font-size: 9px; font-weight: 800; color: #6f746d; letter-spacing: 0.08em; text-transform: uppercase; margin-bottom: 4px;">${escapeHtml(p.label)}</div>
                  <img src="${escapeHtml(p.url)}" alt="${escapeHtml(p.label)}" style="max-width: 100%; max-height: 160px; height: auto; border-radius: 6px; cursor: zoom-in;" onclick="window.open('${escapeHtml(p.url)}', '_blank')">
                </div>
              `).join('')}
            </div>
          ` : ''}

          ${payload.note ? `
            <div style="padding: 10px 18px; background: var(--accent-soft); color: var(--accent-strong); font-size: 12px; font-style: italic;">${escapeHtml(payload.note)}</div>
          ` : ''}

          <!-- Footer -->
          <div style="padding: 10px 18px; background: rgba(245, 241, 232, 0.4); display: flex; justify-content: space-between; align-items: center;">
            <span style="font-size: 11px; color: #6f746d;">${plots.length} plots · ${summary.metric_count || 0} metrics</span>
            <a href="#/jobs/${escapeHtml(compareJobId)}" style="padding: 5px 12px; border-radius: 999px; background: #1c473c; color: #fffdf8; text-decoration: none; font-size: 11px; font-weight: 800;">Open full view →</a>
          </div>
        </div>
      `;
    }).catch(err => {
      wrap.innerHTML = `
        <div style="background: #fffdfa; border: 1px solid #ddd5c7; border-radius: 16px; padding: 18px; max-width: 720px; margin-left: 52px; color: #9b3a2a;">
          Failed to load comparison: ${escapeHtml(String(err))}
        </div>
      `;
    });
    return wrap;
  }

  function showError(msg) {
    errorBanner.textContent = msg;
    errorBanner.style.display = 'block';
    setTimeout(() => { errorBanner.style.display = 'none'; }, 6000);
  }

  async function refreshDrafts() {
    try {
      const drafts = await API.getDesignDrafts(taskId);
      draftCount.textContent = drafts.length;
      draftsList.innerHTML = renderDraftsList(drafts);
      // Toggle showcase link visibility
      const showcaseLink = document.getElementById('design-showcase-link');
      if (showcaseLink) {
        showcaseLink.style.display = drafts.length ? 'inline-flex' : 'none';
      }
      // Wire up delete buttons
      draftsList.querySelectorAll('[data-draft-delete]').forEach(btn => {
        btn.addEventListener('click', async () => {
          const idx = parseInt(btn.getAttribute('data-draft-delete'), 10);
          if (!confirm('Remove this draft?')) return;
          try {
            await API.deleteDesignDraft(taskId, idx);
            await refreshDrafts();
          } catch (e) {
            alert(`Failed: ${e.message}`);
          }
        });
      });
    } catch (e) {
      console.warn('refreshDrafts failed:', e);
    }
  }

  function setupActiveDraftButtons() {
    refreshDrafts();
  }

  function connect() {
    ws = new WebSocket(wsUrl);
    ws.onopen = () => { console.log('Design WS connected'); };
    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      if (data.type === 'session_info') {
        // ignore
      } else if (data.type === 'text_delta') {
        if (!currentAgentMsg) {
          currentAgentMsg = appendMessage('agent', '');
          agentBuffer = '';
        }
        agentBuffer += data.text;
        currentAgentMsg.textContent = agentBuffer;
        chatContainer.scrollTop = chatContainer.scrollHeight;
      } else if (data.type === 'tool_use_start') {
        // Option A: fold any narration + this tool call into the muted
        // "Working" strip instead of breaking the conversation stream.
        foldNarrationIntoProcess();
        addProcessStep('⚙ ' + (data.name || 'tool'));
      } else if (data.type === 'tool_use_end') {
        if (data.ok === false && currentProcessBlock) {
          const body = currentProcessBlock.querySelector('.ap-body');
          if (body && body.lastChild) body.lastChild.textContent += '  ✗';
        }
        if (data.name === 'save_recommendation') {
          refreshDrafts();
        }
      } else if (data.type === 'draft_saved') {
        // Fetch the newly-saved draft and render an inline profile card
        API.getDesignDrafts(taskId).then(drafts => {
          if (drafts && drafts.length) {
            appendMutationProfileCard(drafts[drafts.length - 1]);
          }
        }).catch(err => console.warn('Failed to load new draft:', err));
        refreshDrafts();
      } else if (data.type === 'visualize_residue') {
        // Agent called visualize_residue tool — render an inline 3D viewer
        appendResidueViewerCard(data);
      } else if (data.type === 'comparison_ready') {
        // compare_systems tool finished — show "comparison ready" card
        appendComparisonReadyCard(data);
      } else if (data.type === 'show_comparison_card') {
        // show_comparison_card tool — render an inline diff card
        appendComparisonCard(data);
      } else if (data.type === 'turn_complete') {
        finalizeProcessBlock();
        currentAgentMsg = null;
        agentBuffer = '';
        sendBtn.disabled = false;
        inputEl.disabled = false;
        inputEl.focus();
      } else if (data.type === 'error') {
        finalizeProcessBlock();
        showError(data.message || 'Unknown error');
        sendBtn.disabled = false;
        inputEl.disabled = false;
      }
    };
    ws.onerror = () => { showError('WebSocket error — connection lost'); };
    ws.onclose = () => { console.log('Design WS closed'); };
  }

  function sendMessage(text) {
    if (!text || !text.trim()) return;
    if (!ws || ws.readyState !== WebSocket.OPEN) {
      showError('Not connected. Refresh the page to reconnect.');
      return;
    }
    appendMessage('user', text);
    sendBtn.disabled = true;
    inputEl.disabled = true;
    ws.send(JSON.stringify({ type: 'message', content: text }));
    inputEl.value = '';
  }

  sendBtn.addEventListener('click', () => sendMessage(inputEl.value));
  inputEl.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      sendMessage(inputEl.value);
    }
  });
  document.querySelectorAll('.example-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      const example = btn.getAttribute('data-example');
      if (example) sendMessage(example);
    });
  });

  setupActiveDraftButtons();

  // Pre-render existing drafts as Mutation Profile Cards so the user can
  // see them as soon as they reopen the conversation (rather than only
  // seeing them via the live WebSocket draft_saved event).
  API.getDesignDrafts(taskId).then(drafts => {
    if (!drafts || !drafts.length) return;
    // Insert a small heading message above the cards
    const headerWrap = document.createElement('div');
    headerWrap.style.cssText = 'padding: 12px 16px;';
    headerWrap.innerHTML = `
      <div style="display: flex; align-items: center; gap: 10px; color: #6f746d; font-size: 12px; font-weight: 700; letter-spacing: 0.08em; text-transform: uppercase;">
        <span style="flex: 1; height: 1px; background: #ddd5c7;"></span>
        <span>${drafts.length} previously saved recommendation${drafts.length === 1 ? '' : 's'}</span>
        <span style="flex: 1; height: 1px; background: #ddd5c7;"></span>
      </div>
    `;
    const welcomeMsg = chatContainer.querySelector('.welcome-message');
    if (welcomeMsg) welcomeMsg.remove();
    chatContainer.appendChild(headerWrap);
    drafts.forEach(d => appendMutationProfileCard(d));
  }).catch(err => console.warn('Failed to pre-render drafts:', err));

  // Expose drafts refresh for the multi-agent pipeline panel, and render its
  // recommended sites as rich mutation-profile cards in the conversation when a
  // run completes (the familiar pre-pipeline card form).
  window._designRefreshDrafts = refreshDrafts;
  window._designRenderPipelineCards = async (nAdded) => {
    try {
      const drafts = await API.getDesignDrafts(taskId);
      const recs = (nAdded && nAdded > 0)
        ? drafts.slice(-nAdded)
        : drafts.filter(d => d.source === 'multi_agent');
      if (!recs.length) return;
      const welcome = chatContainer.querySelector('.welcome-message');
      if (welcome) welcome.remove();
      const hdr = document.createElement('div');
      hdr.style.cssText = 'padding: 14px 16px 2px;';
      hdr.innerHTML = `<div style="max-width:720px;margin-left:52px;display:flex;align-items:center;gap:8px;">
        <span style="font-size:9px;font-weight:800;letter-spacing:0.14em;text-transform:uppercase;color:var(--accent-strong);">ImmunoScope · ${recs.length} recommended site${recs.length === 1 ? '' : 's'}</span>
        <span style="flex:1;height:1px;background:var(--line);"></span></div>`;
      chatContainer.appendChild(hdr);
      recs.forEach(d => appendMutationProfileCard(d));
      chatContainer.scrollTop = chatContainer.scrollHeight;
    } catch (e) { console.warn('render pipeline cards failed:', e); }
  };
  // Hydrate the multi-agent pipeline card (left panel) + resume polling if a run
  // is already in progress.
  initMultiAgentPanel(taskId);

  connect();

  // Cleanup on navigation
  window._designCleanup = () => {
    if (ws) { try { ws.close(); } catch(e) {} ws = null; }
    stopMaPolling();
  };
}

// ---------- Showcase: fancy recommendation gallery ----------

// Amino acid properties database (visual fallback when AI doesn't provide one)
const AA_PROPS = {
  A: { name: 'Alanine',       type: 'hydrophobic', polarity: 'nonpolar', charge: 'neutral',  size: 'small',  color: '#94a3b8' },
  R: { name: 'Arginine',      type: 'basic',       polarity: 'polar',    charge: 'positive', size: 'large',  color: '#3b82f6' },
  N: { name: 'Asparagine',    type: 'polar',       polarity: 'polar',    charge: 'neutral',  size: 'medium', color: '#10b981' },
  D: { name: 'Aspartate',     type: 'acidic',      polarity: 'polar',    charge: 'negative', size: 'medium', color: '#ef4444' },
  C: { name: 'Cysteine',      type: 'polar',       polarity: 'polar',    charge: 'neutral',  size: 'small',  color: '#f59e0b' },
  E: { name: 'Glutamate',     type: 'acidic',      polarity: 'polar',    charge: 'negative', size: 'medium', color: '#ef4444' },
  Q: { name: 'Glutamine',     type: 'polar',       polarity: 'polar',    charge: 'neutral',  size: 'medium', color: '#10b981' },
  G: { name: 'Glycine',       type: 'special',     polarity: 'nonpolar', charge: 'neutral',  size: 'small',  color: '#94a3b8' },
  H: { name: 'Histidine',     type: 'basic',       polarity: 'polar',    charge: 'positive', size: 'medium', color: '#3b82f6' },
  I: { name: 'Isoleucine',    type: 'hydrophobic', polarity: 'nonpolar', charge: 'neutral',  size: 'medium', color: '#a78bfa' },
  L: { name: 'Leucine',       type: 'hydrophobic', polarity: 'nonpolar', charge: 'neutral',  size: 'medium', color: '#a78bfa' },
  K: { name: 'Lysine',        type: 'basic',       polarity: 'polar',    charge: 'positive', size: 'large',  color: '#3b82f6' },
  M: { name: 'Methionine',    type: 'hydrophobic', polarity: 'nonpolar', charge: 'neutral',  size: 'medium', color: '#a78bfa' },
  F: { name: 'Phenylalanine', type: 'aromatic',    polarity: 'nonpolar', charge: 'neutral',  size: 'large',  color: '#7c3aed' },
  P: { name: 'Proline',       type: 'special',     polarity: 'nonpolar', charge: 'neutral',  size: 'small',  color: '#f97316' },
  S: { name: 'Serine',        type: 'polar',       polarity: 'polar',    charge: 'neutral',  size: 'small',  color: '#10b981' },
  T: { name: 'Threonine',     type: 'polar',       polarity: 'polar',    charge: 'neutral',  size: 'small',  color: '#10b981' },
  W: { name: 'Tryptophan',    type: 'aromatic',    polarity: 'nonpolar', charge: 'neutral',  size: 'large',  color: '#7c3aed' },
  Y: { name: 'Tyrosine',      type: 'aromatic',    polarity: 'polar',    charge: 'neutral',  size: 'large',  color: '#7c3aed' },
  V: { name: 'Valine',        type: 'hydrophobic', polarity: 'nonpolar', charge: 'neutral',  size: 'medium', color: '#a78bfa' },
};

function getAAProp(code) {
  return AA_PROPS[(code || '').toUpperCase()] || { name: 'Unknown', type: '-', polarity: '-', charge: '-', size: '-', color: '#94a3b8' };
}

function effectIcon(value) {
  const v = String(value || '').toLowerCase();
  // Palette aligned with analysis report: warm earth tones
  if (v.includes('increas') || v.includes('enhanc') || v.includes('improv') || v.includes('gain') || v.includes('强') || v.includes('提升')) return { icon: '↑', color: '#255b4b', bg: '#dbe9e1', label: 'positive' };
  if (v.includes('decreas') || v.includes('reduc') || v.includes('loss') || v.includes('降低') || v.includes('减弱')) return { icon: '↓', color: '#9b3a2a', bg: '#f0ddd5', label: 'negative' };
  if (v.includes('maintain') || v.includes('保持') || v.includes('similar') || v.includes('neutral')) return { icon: '≈', color: '#a36b00', bg: '#f3e8cc', label: 'neutral' };
  return { icon: '·', color: '#6f746d', bg: '#e8e2d4', label: 'unspecified' };
}

// Side-chain chemical features summarized for visual hinting.
// Highlights the most distinctive functional groups for each AA.
const AA_SIDE_CHAIN = {
  A: { features: ['methyl'],                 ring: null,                  charge: null },
  R: { features: ['guanidinium', '+N'],      ring: null,                  charge: '+' },
  N: { features: ['amide', '-CONH2'],        ring: null,                  charge: null },
  D: { features: ['carboxylate', '-COO⁻'],   ring: null,                  charge: '−' },
  C: { features: ['thiol', '-SH'],           ring: null,                  charge: null },
  E: { features: ['carboxylate', '-COO⁻'],   ring: null,                  charge: '−' },
  Q: { features: ['amide', '-CONH2'],        ring: null,                  charge: null },
  G: { features: ['no side chain'],          ring: null,                  charge: null },
  H: { features: ['imidazole', 'N-H'],       ring: 'imidazole (5)',       charge: '±' },
  I: { features: ['branched alkyl'],         ring: null,                  charge: null },
  L: { features: ['branched alkyl'],         ring: null,                  charge: null },
  K: { features: ['lysyl amine', '-NH3⁺'],   ring: null,                  charge: '+' },
  M: { features: ['thioether', '-S-CH3'],    ring: null,                  charge: null },
  F: { features: ['phenyl'],                 ring: 'benzene (6)',         charge: null },
  P: { features: ['cyclic imino'],           ring: 'pyrrolidine (5)',     charge: null },
  S: { features: ['hydroxyl', '-OH'],        ring: null,                  charge: null },
  T: { features: ['hydroxyl', '-OH'],        ring: null,                  charge: null },
  W: { features: ['indole', 'N-H'],          ring: 'indole (5+6)',        charge: null },
  Y: { features: ['phenol', '-OH'],          ring: 'phenol (6)',          charge: null },
  V: { features: ['branched alkyl'],         ring: null,                  charge: null },
};

function getSideChain(code) {
  return AA_SIDE_CHAIN[(code || '').toUpperCase()] || { features: ['unknown'], ring: null, charge: null };
}

// Chemical-property profile for each AA — used to render comparison chips
// in the Mutation Profile Card.
// hydropathy: Kyte-Doolittle scale (-4.5 to +4.5)
// volume: side-chain volume in Å³ (approximate)
// pI: isoelectric point
const AA_PROFILE = {
  A: { type: 'aliphatic', polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: 1.8,  volume: 88.6,  pI: 6.0,  pKa: null },
  R: { type: 'basic',     polarity: 'polar',    charge: 1, hbond: 'donor', hydropathy: -4.5, volume: 173.4, pI: 10.8, pKa: 12.5 },
  N: { type: 'amide',     polarity: 'polar',    charge: 0, hbond: 'both',  hydropathy: -3.5, volume: 114.1, pI: 5.4,  pKa: null },
  D: { type: 'acidic',    polarity: 'polar',    charge: -1, hbond: 'acceptor', hydropathy: -3.5, volume: 111.1, pI: 2.8,  pKa: 3.7 },
  C: { type: 'sulfur',    polarity: 'polar',    charge: 0, hbond: 'donor', hydropathy: 2.5,  volume: 108.5, pI: 5.1,  pKa: 8.3 },
  E: { type: 'acidic',    polarity: 'polar',    charge: -1, hbond: 'acceptor', hydropathy: -3.5, volume: 138.4, pI: 3.2,  pKa: 4.3 },
  Q: { type: 'amide',     polarity: 'polar',    charge: 0, hbond: 'both',  hydropathy: -3.5, volume: 143.8, pI: 5.7,  pKa: null },
  G: { type: 'special',   polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: -0.4, volume: 60.1,  pI: 6.0,  pKa: null },
  H: { type: 'basic',     polarity: 'polar',    charge: 0, hbond: 'both',  hydropathy: -3.2, volume: 153.2, pI: 7.6,  pKa: 6.0 },
  I: { type: 'aliphatic', polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: 4.5,  volume: 166.7, pI: 6.0,  pKa: null },
  L: { type: 'aliphatic', polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: 3.8,  volume: 166.7, pI: 6.0,  pKa: null },
  K: { type: 'basic',     polarity: 'polar',    charge: 1, hbond: 'donor', hydropathy: -3.9, volume: 168.6, pI: 9.7,  pKa: 10.5 },
  M: { type: 'sulfur',    polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: 1.9,  volume: 162.9, pI: 5.7,  pKa: null },
  F: { type: 'aromatic',  polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: 2.8,  volume: 189.9, pI: 5.5,  pKa: null },
  P: { type: 'special',   polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: -1.6, volume: 112.7, pI: 6.3,  pKa: null },
  S: { type: 'hydroxyl',  polarity: 'polar',    charge: 0, hbond: 'both',  hydropathy: -0.8, volume: 89.0,  pI: 5.7,  pKa: null },
  T: { type: 'hydroxyl',  polarity: 'polar',    charge: 0, hbond: 'both',  hydropathy: -0.7, volume: 116.1, pI: 5.6,  pKa: null },
  W: { type: 'aromatic',  polarity: 'nonpolar', charge: 0, hbond: 'donor', hydropathy: -0.9, volume: 227.8, pI: 5.9,  pKa: null },
  Y: { type: 'aromatic',  polarity: 'polar',    charge: 0, hbond: 'both',  hydropathy: -1.3, volume: 193.6, pI: 5.7,  pKa: 10.1 },
  V: { type: 'aliphatic', polarity: 'nonpolar', charge: 0, hbond: '-',    hydropathy: 4.2,  volume: 140.0, pI: 6.0,  pKa: null },
};

function getAAProfile(code) {
  return AA_PROFILE[(code || '').toUpperCase()] || null;
}

// Render a chemical property chip with appropriate color
function aaChip(label, kind = 'default') {
  const colors = {
    aromatic:   { bg: '#ede5f6', color: '#5d4a82', border: '#5d4a8233' },
    aliphatic:  { bg: '#e8e2d4', color: '#6f746d', border: '#6f746d33' },
    polar:      { bg: '#dbe9e1', color: '#255b4b', border: '#255b4b33' },
    nonpolar:   { bg: '#f3e8cc', color: '#a36b00', border: '#a36b0033' },
    acidic:     { bg: '#fde2dc', color: '#9b3a2a', border: '#9b3a2a33' },
    basic:      { bg: '#dde6ef', color: '#3a5f78', border: '#3a5f7833' },
    hydroxyl:   { bg: '#dbe9e1', color: '#255b4b', border: '#255b4b33' },
    amide:      { bg: '#dbe9e1', color: '#255b4b', border: '#255b4b33' },
    sulfur:     { bg: '#fef3c7', color: '#a36b00', border: '#a36b0033' },
    special:    { bg: '#e8e2d4', color: '#6f746d', border: '#6f746d33' },
    plus:       { bg: '#dde6ef', color: '#3a5f78', border: '#3a5f7833' },
    minus:      { bg: '#fde2dc', color: '#9b3a2a', border: '#9b3a2a33' },
    hbond:      { bg: '#ede5f6', color: '#5d4a82', border: '#5d4a8233' },
    default:    { bg: '#fffdfa', color: '#1e2a24', border: '#ddd5c7' },
  };
  const c = colors[kind] || colors.default;
  return `<span style="display: inline-flex; align-items: center; padding: 2px 8px; border-radius: 999px; background: ${c.bg}; color: ${c.color}; border: 1px solid ${c.border}; font-size: 10px; font-weight: 700; letter-spacing: 0.02em; margin: 0 4px 4px 0;">${escapeHtml(label)}</span>`;
}

function aaChipsFor(code) {
  const p = getAAProfile(code);
  if (!p) return '';
  const chips = [];
  chips.push(aaChip(p.type, p.type));
  chips.push(aaChip(p.polarity, p.polarity));
  if (p.charge === 1) chips.push(aaChip('+1 charge', 'plus'));
  else if (p.charge === -1) chips.push(aaChip('−1 charge', 'minus'));
  if (p.hbond && p.hbond !== '-') chips.push(aaChip('H-bond ' + p.hbond, 'hbond'));
  return chips.join('');
}

// Quantitative delta strip — shows numeric changes for the key continuous
// properties (hydropathy, volume) so the user can grasp magnitude of change.
// 3D Viewer modal — uses NGL.js to load the source PDB and highlight the
// mutation residue + its partners.
let _modalEl = null;
function open3DViewer(rec, designTaskId) {
  if (typeof NGL === 'undefined') {
    alert('3D viewer library (NGL) not loaded');
    return;
  }
  // Build the modal once and reuse
  if (_modalEl) _modalEl.remove();
  _modalEl = document.createElement('div');
  _modalEl.style.cssText = `
    position: fixed; inset: 0; z-index: 9999;
    background: rgba(15, 23, 42, 0.65); backdrop-filter: blur(8px);
    display: flex; align-items: center; justify-content: center; padding: 32px;
  `;
  _modalEl.innerHTML = `
    <div style="background: #fffdfa; border-radius: 24px; width: 100%; max-width: 1100px; height: 80vh; display: flex; flex-direction: column; overflow: hidden; box-shadow: 0 30px 80px rgba(0,0,0,0.4);">
      <div style="padding: 16px 22px; border-bottom: 1px solid #ddd5c7; display: flex; justify-content: space-between; align-items: center; background: linear-gradient(180deg, rgba(245, 241, 232, 0.5), transparent);">
        <div>
          <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d;">3D structure viewer</div>
          <div style="font-family: ui-monospace, monospace; font-size: 18px; font-weight: 800; color: #1e2a24; margin-top: 2px;">
            ${escapeHtml(rec.residue || '?')}
            <span style="font-family: 'Noto Sans'; font-size: 13px; font-weight: 500; color: #6f746d; margin-left: 8px;">${escapeHtml(rec.chain || '')} · ${escapeHtml(rec.region || '')}</span>
          </div>
        </div>
        <div style="display: flex; gap: 8px;">
          <select id="ngl-rep-select" style="padding: 6px 10px; border-radius: 8px; border: 1px solid #ddd5c7; background: #fffdfa; font-size: 12px; cursor: pointer;">
            <option value="cartoon">Cartoon</option>
            <option value="surface">Surface</option>
            <option value="ball+stick">Ball &amp; stick</option>
            <option value="ribbon">Ribbon</option>
          </select>
          <button id="ngl-reset" style="padding: 6px 12px; border-radius: 8px; border: 1px solid #ddd5c7; background: #fffdfa; cursor: pointer; font-size: 12px; font-weight: 700;">Reset view</button>
          <button id="ngl-close" style="padding: 6px 14px; border-radius: 8px; border: 0; background: #1c473c; color: #fffdf8; cursor: pointer; font-size: 12px; font-weight: 800;">Close ✕</button>
        </div>
      </div>
      <div id="ngl-stage" style="flex: 1; background: #1e2a24; position: relative;">
        <div id="ngl-loading" style="position: absolute; inset: 0; display: flex; align-items: center; justify-content: center; color: #fffdf8; font-size: 14px;">Loading PDB structure...</div>
      </div>
      <div style="padding: 10px 22px; border-top: 1px solid #ddd5c7; background: #f5f1e8; font-size: 11px; color: #6f746d; display: flex; gap: 18px; flex-wrap: wrap;">
        <span><span style="display: inline-block; width: 12px; height: 12px; border-radius: 50%; background: #9b3a2a; vertical-align: middle; margin-right: 4px;"></span>Mutation residue</span>
        <span><span style="display: inline-block; width: 12px; height: 12px; border-radius: 50%; background: #a36b00; vertical-align: middle; margin-right: 4px;"></span>Contact partners</span>
        <span style="margin-left: auto; font-style: italic;">Drag to rotate · Scroll to zoom · Right-click drag to pan</span>
      </div>
    </div>
  `;
  document.body.appendChild(_modalEl);

  // Close on backdrop click
  _modalEl.addEventListener('click', e => {
    if (e.target === _modalEl) closeModal();
  });
  _modalEl.querySelector('#ngl-close').addEventListener('click', closeModal);
  function closeModal() {
    if (stage) { try { stage.dispose(); } catch(e) {} }
    if (_modalEl) { _modalEl.remove(); _modalEl = null; }
  }

  // Initialize NGL stage
  const stageEl = _modalEl.querySelector('#ngl-stage');
  const stage = new NGL.Stage(stageEl, { backgroundColor: '#1e2a24' });
  window.addEventListener('resize', () => stage.handleResize());

  // PDB endpoint is keyed by the design task id (not the source job id)
  const pdbUrl = `/api/design/jobs/${designTaskId}/structure`;

  Promise.all([
    stage.loadFile(pdbUrl, { ext: 'pdb' }),
    getDesignContextCached(designTaskId),
  ]).then(([component, ctx]) => {
    _modalEl.querySelector('#ngl-loading').remove();

    // Default representation: cartoon
    component.addRepresentation('cartoon', { color: '#cbd5e1', opacity: 0.85 });
    component.autoView();

    // Resolve mutation residue using multi-layer resolver
    const mutLoc = resolveResidueLocation(rec, ctx, false);
    const mutSele = (mutLoc && mutLoc.resid) ? residueSele(mutLoc.resid, mutLoc.chain) : null;
    if (mutLoc) console.log(`[NGL modal] ${rec.residue} → ${mutSele} (via ${mutLoc.source})`);

    if (mutSele) {
      component.addRepresentation('ball+stick', {
        sele: `${mutSele} and not (water or ion)`,
        color: '#9b3a2a',
        radiusScale: 1.6,
      });
      component.addRepresentation('spacefill', {
        sele: mutSele,
        color: '#9b3a2a',
        opacity: 0.35,
      });
    }

    // Highlight partner residues
    const partners = extractContactsFromEvidence(rec);
    const partnerSeles = [];
    partners.forEach(p => {
      const loc = resolvePartnerLocation(p, ctx);
      if (!loc || !loc.resid) return;
      const sele = residueSele(loc.resid, loc.chain);
      partnerSeles.push(sele);
      console.log(`[NGL modal]   partner ${p.residue} → ${sele} (via ${loc.source})`);
      component.addRepresentation('ball+stick', {
        sele: `${sele} and not (water or ion)`,
        color: '#a36b00',
        radiusScale: 1.3,
      });
    });

    // Auto-focus on the mutation residue
    if (mutSele) {
      setTimeout(() => component.autoView(mutSele, 2000), 300);
    }

    // Wire up rep switcher
    const repSelect = _modalEl.querySelector('#ngl-rep-select');
    repSelect.addEventListener('change', e => {
      component.removeAllRepresentations();
      const t = e.target.value;
      component.addRepresentation(t, { color: '#cbd5e1', opacity: t === 'surface' ? 0.5 : 0.85 });
      if (mutSele) {
        component.addRepresentation('ball+stick', { sele: mutSele, color: '#9b3a2a', radiusScale: 1.6 });
        component.addRepresentation('spacefill', { sele: mutSele, color: '#9b3a2a', opacity: 0.35 });
      }
      partnerSeles.forEach(s => component.addRepresentation('ball+stick', { sele: s, color: '#a36b00', radiusScale: 1.3 }));
    });

    _modalEl.querySelector('#ngl-reset').addEventListener('click', () => {
      if (mutSele) component.autoView(mutSele, 1000);
      else component.autoView();
    });
  }).catch(err => {
    const loading = _modalEl.querySelector('#ngl-loading');
    if (loading) {
      loading.innerHTML = `<div style="text-align: center;">
        <p style="font-size: 14px; margin-bottom: 8px;">Failed to load PDB structure</p>
        <p style="font-size: 11px; opacity: 0.7; font-family: ui-monospace, monospace;">${escapeHtml(String(err))}</p>
        <p style="font-size: 11px; opacity: 0.6; margin-top: 12px;">URL: ${escapeHtml(pdbUrl)}</p>
      </div>`;
    }
  });
}

// Parse supporting_evidence md_data lines to extract contact partners.
// Looks for patterns like:
//   "RRCS=8.88 with MET4 (peptide)"
//   "contacts MET4 (peptide), TRP5 (peptide)"
//   "GLN155 (HLA_alpha)"
function extractContactsFromEvidence(rec) {
  const partners = [];
  const ev = rec.supporting_evidence || [];
  for (const item of ev) {
    if (item.type !== 'md_data') continue;
    const text = String(item.data || '');
    // Match "RESIDUE_LABEL (chain_type)" pattern
    // e.g. "MET4 (peptide)", "GLN155 (HLA_alpha)", "TRP5 (peptide)"
    const re = /([A-Z]{3}\d+)\s*\(([a-zA-Z_]+)\)/g;
    let m;
    while ((m = re.exec(text)) !== null) {
      partners.push({ residue: m[1], chain: m[2] });
    }
  }
  // Also check rec.interaction_partners if present (structured form)
  for (const p of (rec.interaction_partners || [])) {
    partners.push({
      residue: p.partner_residue || p.residue || '?',
      chain: p.partner_chain || p.chain || 'peptide',
    });
  }
  return partners;
}

// Counter for unique 3D viewer IDs
let _viewerCounter = 0;

// Cache of design context per task (so 3D viewers can resolve chain IDs).
// Populated lazily by mountInline3DViewer / open3DViewer.
const _designContextCache = {};
async function getDesignContextCached(designTaskId) {
  if (_designContextCache[designTaskId]) return _designContextCache[designTaskId];
  try {
    const ctx = await API.getDesignContext(designTaskId);
    _designContextCache[designTaskId] = ctx;
    return ctx;
  } catch (e) {
    console.warn('Failed to fetch design context for', designTaskId, e);
    return {};
  }
}

// Convert a draft `rec`'s chain field (e.g. "alpha", "TCR_alpha") to the PDB
// chain ID using the chain_mapping from the design context.
function resolveChainId(rawChain, chainMapping) {
  if (!rawChain || !chainMapping) return null;
  const raw = String(rawChain);
  // Greek TCR-domain letters (used in card labels like "TCRβ" / "β") bridge to
  // the word keys the chain_mapping uses ("beta" → "E"). Without this the
  // explicit chain never resolves and the viewer highlights the residue number
  // on EVERY chain.
  const greek = { 'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta' };
  const stripped = raw.replace(/^TCR_?/i, '');  // "TCRβ"→"β", "TCR_alpha"→"alpha"
  const candidates = [
    raw, stripped, stripped.toLowerCase(), raw.toLowerCase(),
    greek[raw], greek[stripped],
  ];
  for (const k of candidates) {
    if (k && chainMapping[k]) return chainMapping[k];
  }
  // Last resort: a bare single-character token is likely already a PDB chain id.
  if (/^[A-Za-z0-9]$/.test(raw)) return raw.toUpperCase();
  return null;
}

// Bullet-proof residue resolver: returns { resid, chain } for a residue
// label like "TYR99", consulting:
//   1. residue_lookup table (PDB-accurate, from RRCS CSV)
//   2. chain_mapping with explicit rec.chain field
//   3. region field (e.g. "CDR3_alpha" → alpha → chain D)
// Returns null if no resolution is possible.
function resolveResidueLocation(rec, ctx, isPartner = false) {
  const ctxLookup = (ctx && ctx.residue_lookup) || {};
  const ctxMapping = (ctx && ctx.chain_mapping) || {};

  const residueLabel = isPartner ? rec.residue : (rec.residue || '');
  if (!residueLabel) return null;

  // ---- Layer 1: direct lookup from RRCS interaction data ----
  if (ctxLookup[residueLabel]) {
    const entry = ctxLookup[residueLabel];
    return { resid: entry.pdb_resid, chain: entry.pdb_chain, source: 'lookup' };
  }

  // Parse number from residue label
  const m = String(residueLabel).match(/\d+/);
  const resNum = m ? parseInt(m[0], 10) : null;
  if (!resNum) return null;

  // ---- Layer 2: explicit chain field ----
  const rawChain = isPartner ? rec.chain : rec.chain;
  let chainId = resolveChainId(rawChain, ctxMapping);
  if (chainId) return { resid: resNum, chain: chainId, source: 'chain_field' };

  // ---- Layer 3: derive chain from region field (e.g. "CDR3_alpha") ----
  if (!isPartner) {
    const region = String(rec.region || '').toLowerCase();
    if (region.includes('alpha')) chainId = ctxMapping['alpha'] || ctxMapping['TCR_alpha'];
    else if (region.includes('beta'))  chainId = ctxMapping['beta']  || ctxMapping['TCR_beta'];
    if (chainId) return { resid: resNum, chain: chainId, source: 'region' };
  }

  // No resolution — return resnum only (Layer 4 fallback)
  return { resid: resNum, chain: null, source: 'fallback' };
}

// Same idea, but for a partner record { residue: "MET4", chain: "peptide" }
function resolvePartnerLocation(partner, ctx) {
  const ctxLookup = (ctx && ctx.residue_lookup) || {};
  const ctxMapping = (ctx && ctx.chain_mapping) || {};
  if (ctxLookup[partner.residue]) {
    const e = ctxLookup[partner.residue];
    return { resid: e.pdb_resid, chain: e.pdb_chain, source: 'lookup' };
  }
  const m = String(partner.residue || '').match(/\d+/);
  const resNum = m ? parseInt(m[0], 10) : null;
  if (!resNum) return null;
  const chainId = resolveChainId(partner.chain, ctxMapping);
  return { resid: resNum, chain: chainId, source: chainId ? 'chain_field' : 'fallback' };
}

// Build NGL selection string for a residue, including chain when available.
// Returns e.g. "99:D" or "99" (fallback when no chain info).
function residueSele(resNum, chainId) {
  if (chainId) return `${resNum}:${chainId}`;
  return `${resNum}`;
}

// Render placeholder for inline 3D viewer. The actual NGL stage is
// initialized in mountInline3DViewer() after the element is in the DOM.
function renderInline3DViewer(rec) {
  const id = `viewer_${++_viewerCounter}`;
  return `
    <div data-inline-viewer="${id}" data-residue="${escapeHtml(rec.residue || '')}"
         style="position: relative; width: 100%; height: 320px; background: #1e2a24; border-radius: 14px; overflow: hidden; box-shadow: inset 0 0 0 1px rgba(255,255,255,0.05);">
      <div class="viewer-loading" style="position: absolute; inset: 0; display: flex; align-items: center; justify-content: center; color: #94a3b8; font-size: 13px; pointer-events: none;">
        <div style="text-align: center;">
          <div style="font-size: 11px; letter-spacing: 0.16em; text-transform: uppercase; opacity: 0.6; margin-bottom: 4px;">Loading PDB structure</div>
          <div class="typing-dots" style="display: flex; gap: 4px; justify-content: center;">
            <span style="width: 6px; height: 6px; border-radius: 50%; background: #94a3b8;"></span>
            <span style="width: 6px; height: 6px; border-radius: 50%; background: #94a3b8;"></span>
            <span style="width: 6px; height: 6px; border-radius: 50%; background: #94a3b8;"></span>
          </div>
        </div>
      </div>
      <!-- Watermark / hint -->
      <div style="position: absolute; bottom: 6px; right: 10px; font-size: 9px; color: rgba(255,255,255,0.35); pointer-events: none; letter-spacing: 0.08em;">drag to rotate · scroll to zoom</div>
    </div>
  `;
}

// Initialize NGL stage for an inline viewer element. Call after the wrap
// is in the DOM. Loads PDB, applies cartoon + highlighted residues, and
// auto-focuses the mutation site.
function mountInline3DViewer(viewerEl, rec, designTaskId) {
  if (typeof NGL === 'undefined') {
    const loading = viewerEl.querySelector('.viewer-loading');
    if (loading) loading.innerHTML = '<div style="color: #b91c1c; font-size: 12px;">NGL library not loaded</div>';
    return;
  }

  // Wait for the element to have a real size before initializing NGL.
  // The viewer element is inside a flex/grid layout inside a scrollable
  // container; on first mount its computed width may still be 0 because
  // the browser hasn't done layout yet. Use rAF + size polling.
  const ensureSizedAndInit = () => {
    const r = viewerEl.getBoundingClientRect();
    if (r.width < 20 || r.height < 20) {
      // Try again next frame (max 30 frames ~= 500ms)
      viewerEl._sizeWaitCount = (viewerEl._sizeWaitCount || 0) + 1;
      if (viewerEl._sizeWaitCount < 30) {
        requestAnimationFrame(ensureSizedAndInit);
        return;
      }
    }
    initStage();
  };

  const loading = viewerEl.querySelector('.viewer-loading');
  let stage = null;

  function initStage() {
  const r0 = viewerEl.getBoundingClientRect();
  console.log(`[NGL inline] init for ${rec.residue || '?'} — container ${r0.width}×${r0.height}`);
  stage = new NGL.Stage(viewerEl, {
    backgroundColor: '#1e2a24',
    quality: 'medium',
    cameraType: 'perspective',
  });
  viewerEl._nglStage = stage;
  // Force a resize after init to make sure canvas matches container
  setTimeout(() => {
    try {
      stage.handleResize();
      const r1 = viewerEl.getBoundingClientRect();
      console.log(`[NGL inline] after resize — ${r1.width}×${r1.height}`);
    } catch (e) { console.warn('NGL resize failed:', e); }
  }, 0);

  // Resolve residue number and partners
  const resMatch = (rec.residue || '').match(/\d+/);
  const resNum = resMatch ? parseInt(resMatch[0], 10) : null;
  const partners = extractContactsFromEvidence(rec);

  // PDB endpoint is keyed by the design task id
  const pdbUrl = `/api/design/jobs/${designTaskId}/structure`;

  // Fetch chain mapping concurrently with PDB load
  const contextPromise = getDesignContextCached(designTaskId);

  Promise.all([
    stage.loadFile(pdbUrl, { ext: 'pdb' }),
    contextPromise,
  ]).then(([component, ctx]) => {
    if (loading) loading.remove();

    // Resolve PDB locations using the multi-layer resolver (lookup → chain → region)
    const mutLoc = resolveResidueLocation(rec, ctx, false);
    const mutSele = (mutLoc && mutLoc.resid) ? residueSele(mutLoc.resid, mutLoc.chain) : null;
    if (mutLoc) console.log(`[NGL] ${rec.residue} → ${mutSele} (resolved via ${mutLoc.source})`);

    // Cartoon backbone (subtle gray-blue)
    component.addRepresentation('cartoon', {
      color: '#7c8a96',
      opacity: 0.75,
      smoothSheet: true,
    });

    // Mutation residue: stick + transparent sphere "halo"
    if (mutSele) {
      component.addRepresentation('ball+stick', {
        sele: `${mutSele} and sidechainAttached`,
        color: '#dc2626',
        radiusScale: 1.5,
        aspectRatio: 1.5,
      });
      // Transparent vdW sphere as visual emphasis
      component.addRepresentation('spacefill', {
        sele: `${mutSele} and sidechainAttached`,
        color: '#dc2626',
        opacity: 0.18,
        radiusScale: 1.0,
      });
      // Label
      component.addRepresentation('label', {
        sele: `${mutSele} and .CA`,
        color: '#fef2f2',
        labelType: 'residue',
        radiusType: 'size',
        radius: 1.2,
        scale: 1.4,
        showBackground: true,
        backgroundColor: '#9b1c1c',
        backgroundOpacity: 0.8,
      });
    }

    // Partner residues: orange sticks
    const partnerSeles = [];
    partners.forEach(p => {
      const loc = resolvePartnerLocation(p, ctx);
      if (!loc || !loc.resid) return;
      const sele = residueSele(loc.resid, loc.chain);
      partnerSeles.push(sele);
      console.log(`[NGL]   partner ${p.residue} → ${sele} (resolved via ${loc.source})`);
      component.addRepresentation('ball+stick', {
        sele: `${sele} and sidechainAttached`,
        color: '#f59e0b',
        radiusScale: 1.1,
        aspectRatio: 1.2,
      });
      component.addRepresentation('label', {
        sele: `${sele} and .CA`,
        color: '#fef3c7',
        labelType: 'residue',
        radiusType: 'size',
        radius: 1.0,
        scale: 1.0,
        showBackground: true,
        backgroundColor: '#92400e',
        backgroundOpacity: 0.7,
      });
    });

    // H-bond representations between mutation residue and partners
    if (mutSele && partnerSeles.length) {
      const allSeles = [mutSele, ...partnerSeles].map(s => `(${s})`).join(' or ');
      component.addRepresentation('contact', {
        sele: `(${allSeles}) and not water`,
        contactType: 'hydrogenBond',
        maxHbondDist: 3.5,
        color: '#fde047',
      });
    }

    // Focus camera on the mutation residue with some context
    setTimeout(() => {
      if (mutSele) {
        try {
          component.autoView(mutSele, 1200);
        } catch (e) {
          component.autoView();
        }
      } else {
        component.autoView();
      }
    }, 200);

    // Handle resize when window changes
    if (!viewerEl._resizeHooked) {
      viewerEl._resizeHooked = true;
      const ro = new ResizeObserver(() => stage.handleResize());
      ro.observe(viewerEl);
      viewerEl._resizeObserver = ro;
    }
  }).catch(err => {
    if (loading) {
      loading.innerHTML = `
        <div style="text-align: center; color: #fca5a5; padding: 12px;">
          <div style="font-size: 13px; font-weight: 700; margin-bottom: 4px;">PDB load failed</div>
          <div style="font-size: 10px; opacity: 0.7; font-family: ui-monospace, monospace;">${escapeHtml(String(err))}</div>
        </div>
      `;
    }
  });
  }  // end initStage

  // Kick off the sized-init dance
  requestAnimationFrame(ensureSizedAndInit);
}

// Render a compact schematic of the TCR-pMHC complex showing where the
// mutation site is and which partner residues it contacts.
function renderComplexSchematic(rec) {
  const W = 680, H = 240;
  const chain = (rec.chain || '').toLowerCase();   // "alpha" | "beta" | "TCR_alpha" etc.
  const region = (rec.region || '').toLowerCase(); // e.g. "cdr3_alpha"
  const isAlpha = chain.includes('alpha') || region.includes('alpha');
  const isBeta = chain.includes('beta') || region.includes('beta');
  const residueLabel = rec.residue || '?';

  // Determine which CDR loop the residue sits on
  let cdrIdx = 3; // default to CDR3 (most likely)
  if (region.includes('cdr1')) cdrIdx = 1;
  else if (region.includes('cdr2')) cdrIdx = 2;
  else if (region.includes('cdr3')) cdrIdx = 3;

  // TCR layout
  const tcrAlphaCx = 200, tcrBetaCx = 480;
  const tcrCy = 70;
  const tcrRx = 130, tcrRy = 44;

  // CDR loop x-positions on the TCR ellipses (3 loops per chain)
  // Order along binding face: CDR1, CDR2, CDR3 (CDR3 in middle)
  const cdrAlphaX = { 1: 130, 2: 200, 3: 270 };
  const cdrBetaX  = { 1: 410, 2: 480, 3: 550 };
  const cdrY = 115; // CDR labels sit at the bottom of TCR ellipse

  // Mutation site coordinates (where the highlighted dot goes)
  const mutX = isAlpha ? cdrAlphaX[cdrIdx] : (isBeta ? cdrBetaX[cdrIdx] : 340);
  const mutY = cdrY + 4;

  // Peptide layout — 9 beads in the groove
  const pepY = 160;
  const pepLeft = 200, pepRight = 480;
  const pepCount = 9;
  const pepBeads = [];
  for (let i = 0; i < pepCount; i++) {
    const x = pepLeft + (pepRight - pepLeft) * (i / (pepCount - 1));
    pepBeads.push({ x, y: pepY, idx: i + 1 });
  }

  // HLA layout — two α-helices (α1 and α2) underneath
  const hlaY = 210;

  // Extract contact partners
  const partners = extractContactsFromEvidence(rec);

  // Try to map partner residues onto schematic positions
  function partnerPos(p) {
    const ch = (p.chain || '').toLowerCase();
    if (ch.includes('pept')) {
      // Try to extract position from residue label (e.g. "MET4" → 4)
      const m = (p.residue || '').match(/\d+/);
      const pos = m ? parseInt(m[0], 10) : Math.floor(pepCount / 2);
      const idx = Math.min(Math.max(pos - 1, 0), pepCount - 1);
      return { x: pepBeads[idx].x, y: pepY, label: p.residue, group: 'peptide' };
    }
    if (ch.includes('hla') || ch.includes('mhc')) {
      // Land somewhere on the HLA helices (alpha1 left, alpha2 right)
      const isA1 = ch.includes('alpha') || ch.includes('a1') || ch.includes('1');
      return { x: isA1 ? 220 : 480, y: hlaY + 8, label: p.residue, group: 'hla' };
    }
    return null;
  }

  // Limit partners to top 3 for visual clarity
  const partnerPositions = partners.slice(0, 3).map(partnerPos).filter(Boolean);

  // SVG
  const palette = {
    tcr: '#dbe9e1',
    tcrStroke: '#255b4b',
    peptide: '#f3e8cc',
    peptideStroke: '#a36b00',
    hla: '#dde6ef',
    hlaStroke: '#3a5f78',
    highlight: '#9b3a2a',
    contact: '#255b4b',
    text: '#1e2a24',
    muted: '#6f746d',
  };

  return `
    <svg viewBox="0 0 ${W} ${H}" width="100%" preserveAspectRatio="xMidYMid meet"
         style="display: block; font-family: 'Noto Sans', sans-serif;">
      <!-- TCR α -->
      <g>
        <ellipse cx="${tcrAlphaCx}" cy="${tcrCy}" rx="${tcrRx}" ry="${tcrRy}"
                 fill="${palette.tcr}" stroke="${palette.tcrStroke}" stroke-width="1.5" opacity="0.85"/>
        <text x="${tcrAlphaCx}" y="${tcrCy - 14}" text-anchor="middle" font-size="11" font-weight="800"
              fill="${palette.tcrStroke}" letter-spacing="0.1em">TCR α</text>
        ${[1, 2, 3].map(i => `
          <g>
            <line x1="${cdrAlphaX[i]}" y1="${tcrCy + tcrRy - 10}" x2="${cdrAlphaX[i]}" y2="${tcrCy + tcrRy + 2}"
                  stroke="${palette.tcrStroke}" stroke-width="1.5"/>
            <text x="${cdrAlphaX[i]}" y="${tcrCy + tcrRy + 14}" text-anchor="middle" font-size="9" font-weight="700"
                  fill="${palette.muted}">CDR${i}</text>
          </g>
        `).join('')}
      </g>
      <!-- TCR β -->
      <g>
        <ellipse cx="${tcrBetaCx}" cy="${tcrCy}" rx="${tcrRx}" ry="${tcrRy}"
                 fill="${palette.tcr}" stroke="${palette.tcrStroke}" stroke-width="1.5" opacity="0.85"/>
        <text x="${tcrBetaCx}" y="${tcrCy - 14}" text-anchor="middle" font-size="11" font-weight="800"
              fill="${palette.tcrStroke}" letter-spacing="0.1em">TCR β</text>
        ${[1, 2, 3].map(i => `
          <g>
            <line x1="${cdrBetaX[i]}" y1="${tcrCy + tcrRy - 10}" x2="${cdrBetaX[i]}" y2="${tcrCy + tcrRy + 2}"
                  stroke="${palette.tcrStroke}" stroke-width="1.5"/>
            <text x="${cdrBetaX[i]}" y="${tcrCy + tcrRy + 14}" text-anchor="middle" font-size="9" font-weight="700"
                  fill="${palette.muted}">CDR${i}</text>
          </g>
        `).join('')}
      </g>

      <!-- Contact lines from mutation site to partners -->
      ${partnerPositions.map(p => `
        <line x1="${mutX}" y1="${mutY}" x2="${p.x}" y2="${p.y}"
              stroke="${palette.contact}" stroke-width="1.5" stroke-dasharray="4 3" opacity="0.6"/>
      `).join('')}

      <!-- Mutation site highlight (drawn AFTER contact lines so it sits on top) -->
      <g>
        <circle cx="${mutX}" cy="${mutY}" r="9" fill="${palette.highlight}" stroke="#fffdfa" stroke-width="2"/>
        <text x="${mutX}" y="${mutY - 14}" text-anchor="middle" font-size="11" font-weight="800"
              fill="${palette.highlight}" font-family="ui-monospace, monospace">${escapeHtml(residueLabel)}</text>
      </g>

      <!-- Peptide groove (background bar) -->
      <rect x="${pepLeft - 20}" y="${pepY - 12}" width="${pepRight - pepLeft + 40}" height="24"
            rx="12" fill="${palette.peptide}" stroke="${palette.peptideStroke}" stroke-width="1.5" opacity="0.85"/>
      <text x="${pepLeft - 50}" y="${pepY + 4}" text-anchor="end" font-size="10" font-weight="800"
            fill="${palette.peptideStroke}" letter-spacing="0.1em">PEPTIDE</text>
      ${pepBeads.map(b => {
        const isPartner = partnerPositions.some(p => p.group === 'peptide' && Math.abs(p.x - b.x) < 3);
        return `
          <circle cx="${b.x}" cy="${b.y}" r="${isPartner ? 7 : 4.5}"
                  fill="${isPartner ? palette.peptideStroke : '#fffdfa'}"
                  stroke="${palette.peptideStroke}" stroke-width="${isPartner ? 2 : 1.4}"/>
          <text x="${b.x}" y="${b.y + 18}" text-anchor="middle" font-size="8" fill="${palette.muted}">${b.idx}</text>
        `;
      }).join('')}

      <!-- HLA helices (two parallel bars) -->
      <g>
        <rect x="60" y="${hlaY}" width="280" height="16" rx="8"
              fill="${palette.hla}" stroke="${palette.hlaStroke}" stroke-width="1.5" opacity="0.85"/>
        <text x="200" y="${hlaY + 11}" text-anchor="middle" font-size="9" font-weight="700"
              fill="${palette.hlaStroke}">HLA α1 helix</text>

        <rect x="340" y="${hlaY}" width="280" height="16" rx="8"
              fill="${palette.hla}" stroke="${palette.hlaStroke}" stroke-width="1.5" opacity="0.85"/>
        <text x="480" y="${hlaY + 11}" text-anchor="middle" font-size="9" font-weight="700"
              fill="${palette.hlaStroke}">HLA α2 helix</text>
      </g>

      <!-- Partner labels -->
      ${partnerPositions.map(p => `
        <g>
          <circle cx="${p.x}" cy="${p.y}" r="7" fill="${p.group === 'peptide' ? palette.peptideStroke : palette.hlaStroke}"
                  stroke="#fffdfa" stroke-width="2"/>
          <text x="${p.x}" y="${p.y + (p.group === 'peptide' ? 30 : 30)}" text-anchor="middle"
                font-size="10" font-weight="800" fill="${p.group === 'peptide' ? palette.peptideStroke : palette.hlaStroke}"
                font-family="ui-monospace, monospace">${escapeHtml(p.label)}</text>
        </g>
      `).join('')}
    </svg>
  `;
}

function renderQuantitativeDeltaStrip(fromCode, toCode) {
  const a = getAAProfile(fromCode);
  const b = getAAProfile(toCode);
  if (!a || !b) return '';

  // Each row: name, fromVal, toVal, unit, axisRange [min, max]
  const rows = [
    { name: 'Hydropathy', from: a.hydropathy, to: b.hydropathy, unit: '', range: [-4.5, 4.5], precision: 1 },
    { name: 'Volume',     from: a.volume,     to: b.volume,     unit: 'Å³', range: [60, 230], precision: 0 },
    { name: 'pI',         from: a.pI,         to: b.pI,         unit: '',  range: [2.5, 11], precision: 1 },
  ];

  return `
    <div style="margin-top: 12px; padding: 12px; background: #fffdfa; border: 1px solid #ddd5c7; border-radius: 12px;">
      <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d; margin-bottom: 10px;">Quantitative shift</div>
      <div style="display: grid; gap: 10px;">
        ${rows.map(row => {
          const pctA = (row.from - row.range[0]) / (row.range[1] - row.range[0]) * 100;
          const pctB = (row.to - row.range[0]) / (row.range[1] - row.range[0]) * 100;
          const left = Math.min(pctA, pctB);
          const width = Math.abs(pctB - pctA);
          const delta = row.to - row.from;
          const deltaSign = delta > 0 ? '+' : '';
          const deltaColor = Math.abs(delta) < (row.range[1] - row.range[0]) * 0.05 ? '#6f746d' : (delta > 0 ? '#255b4b' : '#9b3a2a');
          return `
            <div>
              <div style="display: flex; justify-content: space-between; font-size: 11px; margin-bottom: 4px;">
                <span style="color: #6f746d; font-weight: 600;">${row.name}</span>
                <span style="font-family: ui-monospace, monospace; color: #1e2a24;">
                  <span>${row.from.toFixed(row.precision)}${row.unit}</span>
                  <span style="color: #6f746d;"> → </span>
                  <span style="font-weight: 700;">${row.to.toFixed(row.precision)}${row.unit}</span>
                  <span style="margin-left: 6px; color: ${deltaColor}; font-weight: 700;">(${deltaSign}${delta.toFixed(row.precision)})</span>
                </span>
              </div>
              <div style="position: relative; height: 6px; background: #f3eee2; border-radius: 999px; overflow: visible;">
                <div style="position: absolute; left: ${pctA}%; top: -2px; width: 10px; height: 10px; border-radius: 50%; background: ${a.type === b.type ? '#a36b00' : '#6f746d'}; transform: translateX(-5px); border: 2px solid #fffdfa; box-shadow: 0 0 0 1px #ddd5c7;"></div>
                <div style="position: absolute; left: ${left}%; top: 0; width: ${width}%; height: 100%; background: linear-gradient(90deg, transparent, ${deltaColor}99, transparent); border-radius: 999px;"></div>
                <div style="position: absolute; left: ${pctB}%; top: -3px; width: 12px; height: 12px; border-radius: 50%; background: #255b4b; transform: translateX(-6px); border: 2px solid #fffdfa; box-shadow: 0 2px 4px rgba(37, 91, 75, 0.3);"></div>
              </div>
            </div>
          `;
        }).join('')}
      </div>
    </div>
  `;
}

// SVG side-chain icon — simple chemical topology hint, not atomic-accurate.
// Amino acid SMILES strings (full molecule, in zwitterionic form).
// These are standard SMILES from PubChem / Wikipedia — rendered by
// SmilesDrawer for chemically-correct depiction.
const AA_SMILES = {
  A: 'C[C@@H](N)C(=O)O',                                  // Alanine
  R: 'NC(=N)NCCC[C@@H](N)C(=O)O',                         // Arginine (guanidinium)
  N: 'NC(=O)C[C@@H](N)C(=O)O',                            // Asparagine
  D: 'OC(=O)C[C@@H](N)C(=O)O',                            // Aspartate
  C: 'SC[C@@H](N)C(=O)O',                                 // Cysteine
  E: 'OC(=O)CC[C@@H](N)C(=O)O',                           // Glutamate
  Q: 'NC(=O)CC[C@@H](N)C(=O)O',                           // Glutamine
  G: 'NCC(=O)O',                                          // Glycine
  H: 'O=C(O)[C@@H](N)Cc1cnc[nH]1',                        // Histidine (imidazole)
  I: 'CC[C@H](C)[C@@H](N)C(=O)O',                         // Isoleucine
  L: 'CC(C)C[C@@H](N)C(=O)O',                             // Leucine
  K: 'NCCCC[C@@H](N)C(=O)O',                              // Lysine
  M: 'CSCC[C@@H](N)C(=O)O',                               // Methionine
  F: 'O=C(O)[C@@H](N)Cc1ccccc1',                          // Phenylalanine
  P: 'O=C(O)[C@@H]1CCCN1',                                // Proline (cyclic)
  S: 'OC[C@@H](N)C(=O)O',                                 // Serine
  T: 'C[C@@H](O)[C@@H](N)C(=O)O',                         // Threonine
  W: 'O=C(O)[C@@H](N)Cc1c[nH]c2ccccc12',                  // Tryptophan (indole)
  Y: 'O=C(O)[C@@H](N)Cc1ccc(O)cc1',                       // Tyrosine (phenol)
  V: 'CC(C)[C@@H](N)C(=O)O',                              // Valine
};

// Track SVG containers that need post-render via SmilesDrawer.
// We attach a unique ID, write a placeholder <svg>, and the mount step
// (renderPendingMolecules) walks the DOM and asks SmilesDrawer to draw.
let _smilesCounter = 0;

function renderPendingMolecules(root = document) {
  if (typeof SmilesDrawer === 'undefined') return;
  const targets = root.querySelectorAll('svg[data-smiles]');
  if (!targets.length) return;
  // Ball-and-stick style: colored atom balls + simple bond sticks.
  // Looks like a friendly textbook diagram, not abstract skeleton.
  const drawer = new SmilesDrawer.SvgDrawer({
    width: 160,
    height: 120,
    bondThickness: 2.2,
    bondLength: 22,
    shortBondLength: 0.85,
    bondSpacing: 4.5,
    atomVisualization: 'balls',  // <-- ball-and-stick mode
    isomeric: false,
    debug: false,
    terminalCarbons: true,
    explicitHydrogens: false,
    overlapSensitivity: 0.5,
    overlapResolutionIterations: 5,
    compactDrawing: false,
    fontSizeLarge: 11,
    fontSizeSmall: 6,
    padding: 8,
    experimental: false,
  });
  targets.forEach(svg => {
    const smiles = svg.getAttribute('data-smiles');
    if (!smiles) return;
    if (svg.hasAttribute('data-rendered')) return;
    svg.setAttribute('data-rendered', '1');
    try {
      SmilesDrawer.parse(smiles, tree => {
        try {
          drawer.draw(tree, svg, 'light', false);
        } catch (e) {
          console.warn('SmilesDrawer.draw failed:', smiles, e);
        }
      }, err => {
        console.warn('SmilesDrawer.parse failed:', smiles, err);
      });
    } catch (e) {
      console.warn('SmilesDrawer error for', smiles, e);
    }
  });
}

// Chemically-correct side-chain rendering via SmilesDrawer.
// Returns a placeholder SVG tag with data-smiles attribute; the actual
// chemical structure is drawn by renderPendingMolecules() after the
// element is mounted in the DOM.
// All 20 AAs have Wikipedia-quality structure SVGs locally hosted at
// /static/lib/aa_svg/{code}.svg
const AA_LOCAL_SVG = new Set(['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']);

function aaSideChainSVG(code, color) {
  const c = (code || '').toUpperCase();
  if (!AA_LOCAL_SVG.has(c)) {
    return `<div style="width: 140px; height: 100px; display: flex; align-items: center; justify-content: center; color: #6f746d; font-size: 12px;">${escapeHtml(c)}</div>`;
  }
  return `<img src="/static/lib/aa_svg/${c}.svg" alt="${escapeHtml(c)}" style="width: 100%; max-width: 160px; height: auto; max-height: 120px; object-fit: contain; filter: drop-shadow(0 2px 4px rgba(29,37,32,0.08));" loading="lazy">`;
}

// Compact card for a preserve / do-not-mutate recommendation. The agent
// flagged this residue as load-bearing — there is no substitution to show.
function renderPreserveCard(rec) {
  const wt = (rec.current_aa || '').toUpperCase();
  const wtProp = getAAProp(wt);
  const label = chainRegionLabel(rec) || rec.residue || '?';
  return `
    <div style="background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%); border: 1px solid #ddd5c7; border-radius: 18px; box-shadow: 0 10px 24px rgba(29, 37, 32, 0.06); overflow: hidden; max-width: 720px; margin-left: 52px; font-family: 'Noto Sans', 'PingFang SC', 'Microsoft YaHei', sans-serif; color: #1e2a24;">
      <div style="padding: 12px 18px; background: linear-gradient(180deg, rgba(245, 241, 232, 0.55) 0%, transparent 100%); border-bottom: 1px solid #ddd5c7; display: flex; align-items: center; gap: 10px;">
        <span style="font-size: 14px; font-family: ui-monospace, monospace; font-weight: 800;">${escapeHtml(label)}</span>
        <span style="color: #6f746d; font-size: 11px;">${escapeHtml(rec.residue || '?')}</span>
        <span style="margin-left: auto; padding: 3px 10px; border-radius: 999px; background: #475569; color: #fffdf8; font-size: 9px; font-weight: 800; text-transform: uppercase; letter-spacing: 0.08em;">Preserve</span>
      </div>
      <div style="padding: 16px 18px; display: flex; align-items: center; gap: 14px;">
        <div style="width: 40px; height: 40px; border-radius: 50%; background: ${wtProp.color}; color: #fffdf8; display: flex; align-items: center; justify-content: center; font-family: ui-monospace, monospace; font-size: 17px; font-weight: 800; flex: none;">${escapeHtml(wt || '?')}</div>
        <div>
          <div style="font-size: 13px; font-weight: 700;">Keep ${escapeHtml(wtProp.name || wt)} — do not mutate</div>
          <div class="muted" style="font-size: 11px; margin-top: 2px;">Load-bearing residue; substitution is predicted to hurt binding or specificity.</div>
        </div>
      </div>
      ${rec.rationale ? `<div style="padding: 0 18px 16px; font-size: 12px; line-height: 1.55; color: #44504a;">${escapeHtml(rec.rationale.length > 320 ? rec.rationale.slice(0,320)+'…' : rec.rationale)}</div>` : ''}
    </div>`;
}

// Inline mutation profile card — shown in the chat stream when AI commits
// a recommendation. Designed for at-a-glance "before/after" comprehension.
function renderMutationProfileCard(rec) {
  if (!rec) return '';
  // Preserve / do-not-mutate sites have no real substitution — render a
  // compact "keep this residue" card instead of a meaningless "R → R".
  if (isPreserveDraft(rec)) return renderPreserveCard(rec);
  const currentCode = (rec.current_aa || '').toUpperCase();
  const currentProp = getAAProp(currentCode);
  const suggestions = rec.suggested_mutations || [];
  const firstMut = (suggestions[0] || '').toUpperCase();
  const targetProp = getAAProp(firstMut);
  const currentSC = getSideChain(currentCode);
  const targetSC = getSideChain(firstMut);

  const priorityColors = {
    high: 'linear-gradient(135deg, #1c473c 0%, #2d6f5c 100%)',
    medium: 'linear-gradient(135deg, #a36b00 0%, #d49a4d 100%)',
    low: 'linear-gradient(135deg, #6f746d 0%, #9aa099 100%)',
  };
  const pcBg = priorityColors[rec.priority] || priorityColors.low;

  // Property delta rows
  const propRows = [
    { key: 'type', label: 'Type' },
    { key: 'polarity', label: 'Polarity' },
    { key: 'charge', label: 'Charge' },
    { key: 'size', label: 'Size' },
  ];

  // Extract contacts from supporting_evidence (md_data lines often contain partner info)
  const mdEvidence = (rec.supporting_evidence || []).filter(e => e.type === 'md_data');
  const litEvidence = (rec.supporting_evidence || []).filter(e => e.type === 'literature');

  return `
    <div style="background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%); border: 1px solid #ddd5c7; border-radius: 18px; box-shadow: 0 10px 24px rgba(29, 37, 32, 0.06); overflow: hidden; max-width: 720px; margin-left: 52px; font-family: 'Noto Sans', 'PingFang SC', 'Microsoft YaHei', sans-serif; color: #1e2a24;">
      <!-- Header -->
      <div style="padding: 12px 18px; background: linear-gradient(180deg, rgba(245, 241, 232, 0.55) 0%, transparent 100%); border-bottom: 1px solid #ddd5c7; display: flex; align-items: center; gap: 10px;">
        <span style="font-size: 14px; font-family: ui-monospace, monospace; font-weight: 800;">${escapeHtml(chainRegionLabel(rec) || rec.residue || '?')}</span>
        <span style="color: #6f746d; font-size: 11px;">${escapeHtml(rec.residue || '?')}</span>
        <span style="margin-left: auto; padding: 3px 10px; border-radius: 999px; background: ${pcBg}; color: #fffdf8; font-size: 9px; font-weight: 800; text-transform: uppercase; letter-spacing: 0.08em;">${escapeHtml(rec.priority || 'low')}</span>
        <span style="padding: 3px 10px; border-radius: 999px; background: #fffdfa; border: 1px solid #ddd5c7; color: #1e2a24; font-size: 9px; font-weight: 800; text-transform: uppercase; letter-spacing: 0.08em;">conf: ${escapeHtml(rec.confidence || '-')}</span>
      </div>

      <!-- Before / Arrow / After: ball-and-stick + chips -->
      <div style="padding: 18px 18px; border-bottom: 1px solid #f0ebe0; background: linear-gradient(180deg, rgba(245, 241, 232, 0.30) 0%, transparent 100%);">
        <div style="display: grid; grid-template-columns: 1fr auto 1fr; gap: 14px; align-items: stretch;">
          <!-- BEFORE -->
          <div style="display: flex; flex-direction: column; padding: 14px; border-radius: 14px; background: #fffdfa; border: 1px solid #ddd5c7;">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
              <span style="font-size: 9px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d;">Before</span>
              <div style="width: 32px; height: 32px; border-radius: 50%; background: ${currentProp.color}; color: #fffdf8; display: flex; align-items: center; justify-content: center; font-family: ui-monospace, monospace; font-size: 15px; font-weight: 800;">${escapeHtml(currentCode || '?')}</div>
            </div>
            <div style="display: flex; justify-content: center; align-items: center; min-height: 120px; padding: 4px;">${aaSideChainSVG(currentCode, currentProp.color)}</div>
            <div style="font-size: 13px; font-weight: 700; color: #1e2a24; text-align: center; margin-top: 6px;">${escapeHtml(currentProp.name)}</div>
            <div style="margin-top: 8px; text-align: center; line-height: 1.6;">${aaChipsFor(currentCode)}</div>
          </div>

          <!-- ARROW -->
          <div style="display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 6px;">
            <svg width="44" height="24" viewBox="0 0 44 24" fill="none">
              <path d="M3 12 H36 M30 5 L36 12 L30 19" stroke="#255b4b" stroke-width="2.6" stroke-linecap="round" stroke-linejoin="round"/>
            </svg>
            <span style="font-size: 9px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #6f746d;">mutate</span>
            ${suggestions.length > 1 ? `<span style="font-family: ui-monospace, monospace; font-size: 10px; color: #6f746d;">${suggestions.slice(1).map(escapeHtml).join(' / ')}</span>` : ''}
          </div>

          <!-- AFTER -->
          <div style="display: flex; flex-direction: column; padding: 14px; border-radius: 14px; background: rgba(219, 233, 225, 0.30); border: 1px solid rgba(37, 91, 75, 0.24);">
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
              <span style="font-size: 9px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #255b4b;">After (proposed)</span>
              <div style="width: 32px; height: 32px; border-radius: 50%; background: ${targetProp.color}; color: #fffdf8; display: flex; align-items: center; justify-content: center; font-family: ui-monospace, monospace; font-size: 15px; font-weight: 800;">${escapeHtml(firstMut || '?')}</div>
            </div>
            <div style="display: flex; justify-content: center; align-items: center; min-height: 120px; padding: 4px;">${aaSideChainSVG(firstMut, targetProp.color)}</div>
            <div style="font-size: 13px; font-weight: 700; color: #1e2a24; text-align: center; margin-top: 6px;">${escapeHtml(targetProp.name)}</div>
            <div style="margin-top: 8px; text-align: center; line-height: 1.6;">${aaChipsFor(firstMut)}</div>
          </div>
        </div>

        ${renderQuantitativeDeltaStrip(currentCode, firstMut)}
      </div>

      <!-- Backbone flexibility (φ/ψ from Ramachandran) — pulled from cached context -->
      <div data-backbone-flex-slot="${escapeHtml(rec.residue || '')}" style="display: none;"></div>

      <!-- Property delta table -->
      <div style="padding: 14px 18px; border-bottom: 1px solid #f0ebe0;">
        <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d; margin-bottom: 8px;">Property changes</div>
        <div style="display: grid; grid-template-columns: 80px 1fr 24px 1fr 50px; row-gap: 4px; column-gap: 8px; font-size: 12px; align-items: center;">
          ${propRows.map(row => {
            const same = currentProp[row.key] === targetProp[row.key];
            const marker = same
              ? '<span style="color: #255b4b; font-weight: 800;">✓ same</span>'
              : '<span style="color: #a36b00; font-weight: 800;">△ shift</span>';
            const cellBg = same ? '#fffdfa' : '#f3e8cc';
            return `
              <div style="color: #6f746d; font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em;">${row.label}</div>
              <div style="text-align: center; padding: 4px 8px; background: ${cellBg}; border-radius: 6px; color: #1e2a24;">${escapeHtml(currentProp[row.key])}</div>
              <div style="text-align: center; color: #6f746d;">→</div>
              <div style="text-align: center; padding: 4px 8px; background: ${cellBg}; border-radius: 6px; color: #1e2a24; ${same ? '' : 'font-weight: 700;'}">${escapeHtml(targetProp[row.key])}</div>
              <div style="font-size: 10px; text-align: right;">${marker}</div>
            `;
          }).join('')}
        </div>
      </div>

      <!-- Inline 3D viewer (PyMOL-style structural view) -->
      <div style="padding: 14px 18px; border-bottom: 1px solid #f0ebe0;">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
          <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d;">Structural context</div>
          <button class="btn-view-3d" data-rec="${escapeHtml(rec.residue || '?')}" style="padding: 4px 10px; border-radius: 999px; background: #fffdfa; border: 1px solid #ddd5c7; color: #255b4b; font-size: 10px; font-weight: 800; cursor: pointer; font-family: 'Noto Sans', sans-serif; letter-spacing: 0.06em;">Expand →</button>
        </div>
        ${renderInline3DViewer(rec)}
        <div style="display: flex; justify-content: center; gap: 18px; margin-top: 8px; font-size: 10px; color: #6f746d;">
          <span style="display: inline-flex; align-items: center; gap: 5px;"><span style="display: inline-block; width: 10px; height: 10px; border-radius: 50%; background: #dc2626;"></span>Mutation residue</span>
          <span style="display: inline-flex; align-items: center; gap: 5px;"><span style="display: inline-block; width: 10px; height: 10px; border-radius: 50%; background: #f59e0b;"></span>Contact partner</span>
          <span style="display: inline-flex; align-items: center; gap: 5px;"><span style="display: inline-block; width: 14px; border-top: 1.5px dashed #fde047;"></span>H-bond</span>
        </div>
      </div>

      <!-- Contact environment summary -->
      ${mdEvidence.length || litEvidence.length ? `
        <div style="padding: 14px 18px; border-bottom: 1px solid #f0ebe0;">
          <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d; margin-bottom: 8px;">Contact environment</div>
          ${mdEvidence.length ? `
            <div style="display: flex; flex-direction: column; gap: 4px; font-size: 12px;">
              ${mdEvidence.slice(0, 3).map(e => `
                <div style="display: flex; gap: 8px; align-items: center;">
                  <span style="width: 6px; height: 6px; border-radius: 50%; background: #255b4b; flex-shrink: 0;"></span>
                  <span style="font-family: ui-monospace, monospace; color: #1e2a24;">${escapeHtml(e.data || '-')}</span>
                </div>
              `).join('')}
            </div>
          ` : ''}
          ${litEvidence.length ? `
            <div style="margin-top: 8px; display: flex; flex-wrap: wrap; gap: 4px;">
              ${litEvidence.slice(0, 3).map(e => e.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${escapeHtml(e.pmid)}/" target="_blank" rel="noopener" style="padding: 2px 8px; background: #a36b00; color: #fffdf8; border-radius: 6px; font-size: 10px; font-weight: 800; text-decoration: none;">PMID ${escapeHtml(e.pmid)}</a>` : '').join('')}
            </div>
          ` : ''}
        </div>
      ` : ''}

      <!-- Risk -->
      ${rec.risks && rec.risks.length ? `
        <div style="padding: 12px 18px; background: linear-gradient(180deg, #f9efe9, #fffdfa); display: flex; gap: 10px; align-items: flex-start;">
          <span style="width: 18px; height: 18px; border-radius: 50%; background: rgba(155, 58, 42, 0.15); color: #9b3a2a; text-align: center; line-height: 18px; font-size: 11px; font-weight: 800; flex-shrink: 0; margin-top: 1px;">!</span>
          <div style="font-size: 12px; color: #1e2a24; line-height: 1.5;">
            <strong style="color: #9b3a2a;">Risk:</strong>
            ${escapeHtml(rec.risks[0])}
          </div>
        </div>
      ` : ''}

      <!-- Footer link -->
      <div style="padding: 8px 18px; background: #fbf8f2; display: flex; justify-content: space-between; align-items: center; font-size: 11px; color: #6f746d;">
        <span>Saved to draft panel · #${(rec._draft_index != null ? rec._draft_index + 1 : '?')}</span>
        <span style="font-style: italic;">Click "View Showcase" in the right panel for full detail.</span>
      </div>
    </div>
  `;
}

async function renderDesignShowcase(taskId) {
  let task, drafts, context;
  try {
    const data = await API.getDesignJob(taskId);
    task = data.task;
    drafts = data.drafts || [];
    context = await API.getDesignContext(taskId);
  } catch (e) {
    return `<div class="panel"><h2>Error loading showcase</h2><p>${escapeHtml(e.message)}</p><a class="button" href="#/copilot">Back</a></div>`;
  }

  if (!drafts.length) {
    return `
      <div class="panel" style="text-align: center; padding: 60px 20px;">
        <h2>No recommendations yet</h2>
        <p class="muted">Go back to the conversation and let the AI investigate hotspots and save recommendations.</p>
        <a class="button" href="#/copilot/${escapeHtml(taskId)}">← Back to conversation</a>
      </div>
    `;
  }

  // Summary stats
  const priorities = { high: 0, medium: 0, low: 0 };
  const confidences = { high: 0, medium: 0, low: 0 };
  const regionCounts = {};
  drafts.forEach(d => {
    if (d.priority) priorities[d.priority] = (priorities[d.priority] || 0) + 1;
    if (d.confidence) confidences[d.confidence] = (confidences[d.confidence] || 0) + 1;
    const r = d.region || 'Unknown';
    regionCounts[r] = (regionCounts[r] || 0) + 1;
  });

  // Themed wrapper - all inner colors aligned with analysis report:
  // bg #f5f1e8, paper #fffdfa, ink #1e2a24, accent #255b4b, warn #a36b00
  return `
    <div id="showcase-root" style="font-family: 'Noto Sans', 'PingFang SC', 'Microsoft YaHei', sans-serif; color: #1e2a24; max-width: 1200px; margin: 0 auto; padding: 0 4px;">
      ${renderShowcaseHeader(task, context, drafts, priorities, confidences, regionCounts)}
      <div id="showcase-cards" style="display: flex; flex-direction: column; gap: 28px; margin-top: 28px;">
        ${drafts.map((d, idx) => renderShowcaseCard(d, idx, drafts.length)).join('')}
      </div>
      <div style="display: flex; justify-content: center; gap: 12px; margin-top: 40px; padding-bottom: 40px;">
        <button id="showcase-export-md" style="border: 0; border-radius: 999px; padding: 11px 22px; background: linear-gradient(135deg, #1c473c 0%, #2d6f5c 100%); color: #f7f5ef; font-weight: 700; cursor: pointer; box-shadow: 0 8px 18px rgba(28, 71, 60, 0.18); font-family: inherit;">Export as Markdown</button>
        <button id="showcase-export-json" style="border: 1px solid #ddd5c7; border-radius: 999px; padding: 11px 22px; background: #fffdfa; color: #255b4b; font-weight: 700; cursor: pointer; font-family: inherit;">Export as JSON</button>
        <a href="#/copilot/${escapeHtml(taskId)}" style="display: inline-flex; align-items: center; padding: 11px 22px; border: 1px solid #ddd5c7; border-radius: 999px; background: #fffdfa; color: #6f746d; text-decoration: none; font-weight: 700; font-family: inherit;">← Back to conversation</a>
      </div>
    </div>
  `;
}

function renderShowcaseHeader(task, context, drafts, priorities, confidences, regionCounts) {
  const topRegion = Object.entries(regionCounts).sort((a, b) => b[1] - a[1])[0];
  return `
    <div style="background: linear-gradient(135deg, #1c473c 0%, #2d6f5c 100%); border-radius: 24px; padding: 36px 40px; color: #f7f5ef; position: relative; overflow: hidden; box-shadow: 0 18px 44px rgba(28, 71, 60, 0.18);">
      <!-- Subtle warm decorative overlays (gold + ivory) -->
      <div style="position: absolute; top: -60px; right: -50px; width: 240px; height: 240px; border-radius: 50%; background: radial-gradient(circle at center, rgba(255, 209, 102, 0.18) 0%, transparent 70%); pointer-events: none;"></div>
      <div style="position: absolute; bottom: -90px; left: -80px; width: 280px; height: 280px; border-radius: 50%; background: rgba(247, 245, 239, 0.05); pointer-events: none;"></div>

      <div style="position: relative;">
        <div style="display: flex; justify-content: space-between; align-items: flex-start; gap: 24px; flex-wrap: wrap;">
          <div style="flex: 1; min-width: 240px;">
            <span style="display: inline-block; padding: 4px 12px; border-radius: 999px; background: rgba(255, 209, 102, 0.18); color: #ffd166; font-size: 11px; font-weight: 700; letter-spacing: 0.16em; text-transform: uppercase;">Mutation Design Recommendations</span>
            <h1 style="margin: 14px 0 4px; font-size: 32px; font-weight: 800; letter-spacing: -0.02em; color: #fffdf8;">${escapeHtml(task.name)}</h1>
            <p style="margin: 0; opacity: 0.85; font-size: 13px;">Generated from <strong style="font-family: ui-monospace, SFMono-Regular, Menlo, monospace;">${escapeHtml(task.source_job_name)}</strong></p>
          </div>
          <div style="display: flex; align-items: center; gap: 24px;">
            <div style="text-align: center;">
              <div style="font-size: 38px; font-weight: 800; line-height: 1; color: #fffdf8;">${drafts.length}</div>
              <div style="font-size: 10px; letter-spacing: 0.1em; text-transform: uppercase; opacity: 0.8; margin-top: 4px;">Recommendations</div>
            </div>
            <div style="width: 1px; height: 50px; background: rgba(247, 245, 239, 0.25);"></div>
            <div style="text-align: center;">
              <div style="font-size: 38px; font-weight: 800; line-height: 1; color: #ffd166;">${priorities.high}</div>
              <div style="font-size: 10px; letter-spacing: 0.1em; text-transform: uppercase; opacity: 0.8; margin-top: 4px;">High priority</div>
            </div>
          </div>
        </div>

        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; margin-top: 28px;">
          <div style="background: rgba(247, 245, 239, 0.10); backdrop-filter: blur(10px); border-radius: 18px; padding: 14px 18px; border: 1px solid rgba(247, 245, 239, 0.12);">
            <div style="font-size: 10px; letter-spacing: 0.12em; text-transform: uppercase; opacity: 0.85; font-weight: 700;">System</div>
            <div style="font-size: 14px; font-weight: 600; margin-top: 6px; color: #fffdf8;">${escapeHtml(context.peptide || context.system_id || '-')}</div>
            ${context.hla ? `<div style="font-size: 11px; opacity: 0.85; margin-top: 2px;">${escapeHtml(context.hla)}</div>` : ''}
          </div>
          <div style="background: rgba(247, 245, 239, 0.10); backdrop-filter: blur(10px); border-radius: 18px; padding: 14px 18px; border: 1px solid rgba(247, 245, 239, 0.12);">
            <div style="font-size: 10px; letter-spacing: 0.12em; text-transform: uppercase; opacity: 0.85; font-weight: 700;">Simulation</div>
            <div style="font-size: 14px; font-weight: 600; margin-top: 6px; color: #fffdf8;">${context.duration_ns || '-'} ns · ${context.n_frames || '-'} frames</div>
            <div style="font-size: 11px; opacity: 0.85; margin-top: 2px;">RMSD ${context.rmsd || '-'} nm · ${escapeHtml(context.convergence || '-')}</div>
          </div>
          <div style="background: rgba(247, 245, 239, 0.10); backdrop-filter: blur(10px); border-radius: 18px; padding: 14px 18px; border: 1px solid rgba(247, 245, 239, 0.12);">
            <div style="font-size: 10px; letter-spacing: 0.12em; text-transform: uppercase; opacity: 0.85; font-weight: 700;">Confidence mix</div>
            <div style="display: flex; gap: 4px; margin-top: 8px; height: 8px; border-radius: 999px; overflow: hidden;">
              ${confidences.high ? `<div style="flex: ${confidences.high}; background: #7ec2ad;" title="High: ${confidences.high}"></div>` : ''}
              ${confidences.medium ? `<div style="flex: ${confidences.medium}; background: #ffd166;" title="Medium: ${confidences.medium}"></div>` : ''}
              ${confidences.low ? `<div style="flex: ${confidences.low}; background: #b8b3a4;" title="Low: ${confidences.low}"></div>` : ''}
            </div>
            <div style="font-size: 11px; opacity: 0.85; margin-top: 6px;">H:${confidences.high || 0} · M:${confidences.medium || 0} · L:${confidences.low || 0}</div>
          </div>
          <div style="background: rgba(247, 245, 239, 0.10); backdrop-filter: blur(10px); border-radius: 18px; padding: 14px 18px; border: 1px solid rgba(247, 245, 239, 0.12);">
            <div style="font-size: 10px; letter-spacing: 0.12em; text-transform: uppercase; opacity: 0.85; font-weight: 700;">Top region</div>
            <div style="font-size: 14px; font-weight: 600; margin-top: 6px; color: #fffdf8;">${escapeHtml(topRegion ? topRegion[0] : '-')}</div>
            <div style="font-size: 11px; opacity: 0.85; margin-top: 2px;">${topRegion ? `${topRegion[1]} recommendation${topRegion[1] === 1 ? '' : 's'}` : '-'}</div>
          </div>
        </div>
      </div>
    </div>
  `;
}

function renderShowcaseCard(rec, idx, total) {
  const priorityColors = {
    high: { bg: 'linear-gradient(135deg, #1c473c 0%, #2d6f5c 100%)', label: '#1c473c', text: '#f7f5ef' },
    medium: { bg: 'linear-gradient(135deg, #a36b00 0%, #d49a4d 100%)', label: '#a36b00', text: '#fffdf8' },
    low: { bg: 'linear-gradient(135deg, #6f746d 0%, #9aa099 100%)', label: '#6f746d', text: '#fffdf8' },
  };
  const pc = priorityColors[rec.priority] || priorityColors.low;
  const currentProp = getAAProp(rec.current_aa);

  return `
    <article style="background: linear-gradient(180deg, #ffffff 0%, #fbfdfc 100%); border-radius: 24px; box-shadow: 0 12px 30px rgba(29, 37, 32, 0.08); overflow: hidden; border: 1px solid #ddd5c7;">
      <!-- Card hero -->
      <div style="padding: 28px 32px; display: flex; justify-content: space-between; align-items: flex-start; gap: 20px; border-bottom: 1px solid #ddd5c7; background: linear-gradient(180deg, rgba(245, 241, 232, 0.55) 0%, transparent 100%);">
        <div style="flex: 1; min-width: 0;">
          <div style="display: flex; align-items: center; gap: 10px; margin-bottom: 8px; flex-wrap: wrap;">
            <span style="font-size: 11px; font-weight: 800; color: #6f746d; letter-spacing: 0.14em; text-transform: uppercase;">#${idx + 1} of ${total}</span>
            <span style="padding: 4px 12px; border-radius: 999px; background: ${pc.bg}; color: ${pc.text}; font-size: 10px; font-weight: 800; text-transform: uppercase; letter-spacing: 0.08em;">${escapeHtml(rec.priority || 'low')} priority</span>
            <span style="padding: 4px 12px; border-radius: 999px; background: #fffdfa; border: 1px solid #ddd5c7; color: #1e2a24; font-size: 10px; font-weight: 800; text-transform: uppercase; letter-spacing: 0.08em;">conf: ${escapeHtml(rec.confidence || '-')}</span>
          </div>
          <h2 style="margin: 0; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 36px; font-weight: 800; letter-spacing: -0.02em; line-height: 1; color: #1e2a24;">${escapeHtml(rec.residue || '?')}</h2>
          <p style="margin: 8px 0 0; color: #6f746d; font-size: 14px;">${escapeHtml(chainRegionLabel(rec) || ((rec.chain || '?') + ' chain'))}</p>
        </div>

        <!-- Mutation transformation -->
        <div style="display: flex; align-items: center; gap: 12px; flex-shrink: 0;">
          ${renderAABubble(rec.current_aa, currentProp, true)}
          <div style="display: flex; flex-direction: column; align-items: center; color: #6f746d;">
            <svg width="50" height="20" viewBox="0 0 50 20" fill="none" xmlns="http://www.w3.org/2000/svg">
              <path d="M5 10 H40 M35 5 L40 10 L35 15" stroke="${pc.label}" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
            </svg>
            <span style="font-size: 9px; font-weight: 700; letter-spacing: 0.1em; text-transform: uppercase; margin-top: 4px; color: #6f746d;">mutate to</span>
          </div>
          <div style="display: flex; gap: 8px;">
            ${(rec.suggested_mutations || []).map(m => renderAABubble(m, getAAProp(m), false)).join('')}
          </div>
        </div>
      </div>

      <!-- Body sections -->
      <div style="padding: 28px 32px;">
        <!-- Rationale -->
        ${rec.rationale ? `
          <section style="margin-bottom: 28px;">
            <h3 style="margin: 0 0 12px; font-size: 12px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #255b4b;">Design Rationale</h3>
            <blockquote style="margin: 0; padding: 16px 20px; border-left: 3px solid #255b4b; background: linear-gradient(90deg, rgba(37, 91, 75, 0.05), transparent); border-radius: 0 12px 12px 0; font-size: 14px; line-height: 1.7; color: #1e2a24;">${escapeHtml(rec.rationale)}</blockquote>
          </section>
        ` : ''}

        <!-- Expected effects grid -->
        ${rec.expected_effects && Object.keys(rec.expected_effects).length ? `
          <section style="margin-bottom: 28px;">
            <h3 style="margin: 0 0 14px; font-size: 12px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #255b4b;">Expected Effects</h3>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px;">
              ${Object.entries(rec.expected_effects).map(([key, val]) => {
                const e = effectIcon(val);
                return `
                  <div style="background: #fffdfa; border: 1px solid #ddd5c7; border-radius: 14px; padding: 14px; display: flex; gap: 12px; align-items: flex-start;">
                    <div style="width: 36px; height: 36px; border-radius: 50%; background: ${e.bg}; color: ${e.color}; display: flex; align-items: center; justify-content: center; font-size: 18px; font-weight: 800; flex-shrink: 0;">${e.icon}</div>
                    <div style="flex: 1; min-width: 0;">
                      <div style="font-size: 11px; font-weight: 800; color: #6f746d; text-transform: uppercase; letter-spacing: 0.06em;">${escapeHtml(key)}</div>
                      <div style="font-size: 13px; margin-top: 3px; line-height: 1.4; color: #1e2a24;">${escapeHtml(val)}</div>
                    </div>
                  </div>
                `;
              }).join('')}
            </div>
          </section>
        ` : ''}

        <!-- Evidence: 2-column -->
        <section style="margin-bottom: 28px;">
          <h3 style="margin: 0 0 14px; font-size: 12px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #255b4b;">Supporting Evidence</h3>
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 14px;">
            ${renderEvidenceColumn(rec.supporting_evidence || [], 'md_data', 'MD Data', '#255b4b')}
            ${renderEvidenceColumn(rec.supporting_evidence || [], 'literature', 'Literature', '#a36b00')}
          </div>
        </section>

        <!-- Risks + Validation: 2-column -->
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 18px;">
          <!-- Risks -->
          ${rec.risks && rec.risks.length ? `
            <section style="background: linear-gradient(180deg, #f9efe9, #fffdfa); border: 1px solid #e8c8b8; border-radius: 16px; padding: 18px;">
              <h3 style="margin: 0 0 12px; font-size: 11px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #9b3a2a; display: flex; align-items: center; gap: 8px;">
                <span style="display: inline-block; width: 22px; height: 22px; border-radius: 50%; background: rgba(155, 58, 42, 0.12); color: #9b3a2a; text-align: center; line-height: 22px; font-size: 13px; font-weight: 800;">!</span>
                Risks &amp; Caveats
              </h3>
              <ul style="margin: 0; padding-left: 20px; font-size: 13px; line-height: 1.6; color: #1e2a24;">
                ${rec.risks.map(r => `<li style="margin-bottom: 4px;">${escapeHtml(r)}</li>`).join('')}
              </ul>
            </section>
          ` : '<div></div>'}

          <!-- Validation -->
          ${rec.validation_experiments && rec.validation_experiments.length ? `
            <section style="background: linear-gradient(180deg, #f0eee5, #fffdfa); border: 1px solid #c9c2af; border-radius: 16px; padding: 18px;">
              <h3 style="margin: 0 0 12px; font-size: 11px; font-weight: 800; letter-spacing: 0.12em; text-transform: uppercase; color: #3a5f78; display: flex; align-items: center; gap: 8px;">
                <span style="display: inline-block; width: 22px; height: 22px; border-radius: 50%; background: rgba(58, 95, 120, 0.12); color: #3a5f78; text-align: center; line-height: 22px; font-size: 13px; font-weight: 800;">✓</span>
                Suggested Validation
              </h3>
              <ul style="margin: 0; padding-left: 20px; font-size: 13px; line-height: 1.6; color: #1e2a24;">
                ${rec.validation_experiments.map(v => `<li style="margin-bottom: 4px;">${escapeHtml(v)}</li>`).join('')}
              </ul>
            </section>
          ` : '<div></div>'}
        </div>
      </div>

      <!-- Footer with AA property comparison -->
      ${renderAAComparisonStrip(rec)}
    </article>
  `;
}

function renderAABubble(code, prop, isCurrent) {
  const size = isCurrent ? 72 : 60;
  const fontSize = isCurrent ? 28 : 22;
  return `
    <div style="position: relative; width: ${size}px; height: ${size}px; flex-shrink: 0;" title="${escapeHtml(prop.name)} — ${escapeHtml(prop.type)}, ${escapeHtml(prop.polarity)}">
      <div style="width: 100%; height: 100%; border-radius: 50%; background: ${prop.color}; display: flex; align-items: center; justify-content: center; color: #fffdf8; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: ${fontSize}px; font-weight: 800; box-shadow: 0 8px 20px ${prop.color}33; border: 2px solid #fffdfa;">
        ${escapeHtml(code || '?')}
      </div>
      ${isCurrent ? `<div style="position: absolute; bottom: -14px; left: 50%; transform: translateX(-50%); font-size: 9px; font-weight: 800; letter-spacing: 0.1em; text-transform: uppercase; color: #6f746d; white-space: nowrap;">current</div>` : ''}
    </div>
  `;
}

function renderEvidenceColumn(evidence, type, title, accentColor) {
  const items = evidence.filter(e => e.type === type);
  if (!items.length) {
    return `
      <div style="background: #f5f1e8; border: 1px dashed #ddd5c7; border-radius: 14px; padding: 18px; text-align: center;">
        <div style="font-size: 11px; font-weight: 800; letter-spacing: 0.1em; text-transform: uppercase; color: #6f746d;">${escapeHtml(title)}</div>
        <p style="margin: 8px 0 0; font-size: 12px; color: #8a8f88;">No ${type.replace('_', ' ')} evidence provided</p>
      </div>
    `;
  }
  return `
    <div style="background: #fffdfa; border: 1px solid #ddd5c7; border-radius: 14px; padding: 16px;">
      <div style="display: flex; align-items: center; gap: 8px; margin-bottom: 12px;">
        <span style="width: 8px; height: 8px; border-radius: 50%; background: ${accentColor};"></span>
        <span style="font-size: 11px; font-weight: 800; letter-spacing: 0.1em; text-transform: uppercase; color: ${accentColor};">${escapeHtml(title)}</span>
        <span style="margin-left: auto; padding: 2px 8px; border-radius: 999px; background: ${accentColor}1a; color: ${accentColor}; font-size: 10px; font-weight: 700;">${items.length}</span>
      </div>
      <div style="display: flex; flex-direction: column; gap: 10px;">
        ${items.map(item => {
          if (type === 'md_data') {
            return `<div style="font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; padding: 10px 12px; background: #f5f1e8; border-radius: 8px; color: #1e2a24;">${escapeHtml(item.data || '-')}</div>`;
          } else {
            const pmid = item.pmid || '';
            const title2 = item.title || '-';
            const relevance = item.relevance || '';
            return `
              <div style="padding: 10px 12px; background: #f5f1e8; border-left: 2px solid ${accentColor}; border-radius: 0 8px 8px 0;">
                <div style="display: flex; gap: 8px; align-items: flex-start;">
                  ${pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${escapeHtml(pmid)}/" target="_blank" rel="noopener" style="flex-shrink: 0; padding: 2px 8px; background: ${accentColor}; color: #fffdf8; border-radius: 6px; font-size: 10px; font-weight: 800; text-decoration: none; letter-spacing: 0.04em;">PMID ${escapeHtml(pmid)}</a>` : ''}
                  <div style="flex: 1; font-size: 12px; line-height: 1.5; color: #1e2a24;">
                    <div style="font-weight: 600;">${escapeHtml(title2)}</div>
                    ${relevance ? `<div style="margin-top: 4px; font-size: 11px; color: #6f746d;">${escapeHtml(relevance)}</div>` : ''}
                  </div>
                </div>
              </div>
            `;
          }
        }).join('')}
      </div>
    </div>
  `;
}

function renderAAComparisonStrip(rec) {
  const current = getAAProp(rec.current_aa);
  const suggested = (rec.suggested_mutations || []).map(m => ({ code: m, prop: getAAProp(m) }));
  if (!suggested.length) return '';

  const propRows = [
    { key: 'type',     label: 'Type' },
    { key: 'polarity', label: 'Polarity' },
    { key: 'charge',   label: 'Charge' },
    { key: 'size',     label: 'Size' },
  ];

  return `
    <div style="border-top: 1px solid #ddd5c7; background: #f5f1e8; padding: 18px 32px;">
      <div style="font-size: 10px; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; color: #6f746d; margin-bottom: 12px;">Amino Acid Property Comparison</div>
      <div style="display: grid; grid-template-columns: 120px repeat(${suggested.length + 1}, minmax(0, 1fr)); gap: 8px; font-size: 12px;">
        <div></div>
        <div style="text-align: center;">
          <div style="display: inline-block; padding: 4px 12px; border-radius: 999px; background: ${current.color}22; color: ${current.color}; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-weight: 800;">${escapeHtml(rec.current_aa || '?')} (current)</div>
        </div>
        ${suggested.map(s => `
          <div style="text-align: center;">
            <div style="display: inline-block; padding: 4px 12px; border-radius: 999px; background: ${s.prop.color}22; color: ${s.prop.color}; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-weight: 800;">${escapeHtml(s.code)}</div>
          </div>
        `).join('')}

        ${propRows.map(row => `
          <div style="color: #6f746d; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em; font-size: 10px; align-self: center;">${row.label}</div>
          <div style="text-align: center; padding: 6px; background: #fffdfa; border-radius: 8px; color: #1e2a24;">${escapeHtml(current[row.key])}</div>
          ${suggested.map(s => {
            const same = s.prop[row.key] === current[row.key];
            return `<div style="text-align: center; padding: 6px; background: ${same ? '#fffdfa' : '#f3e8cc'}; border-radius: 8px; color: #1e2a24; ${same ? '' : 'font-weight: 700;'}">${escapeHtml(s.prop[row.key])}</div>`;
          }).join('')}
        `).join('')}
      </div>
    </div>
  `;
}

function mountDesignShowcase(taskId) {
  // Export as Markdown
  const mdBtn = document.getElementById('showcase-export-md');
  if (mdBtn) {
    mdBtn.addEventListener('click', async () => {
      const drafts = await API.getDesignDrafts(taskId);
      const data = await API.getDesignJob(taskId);
      const task = data.task;
      const md = generateMarkdownReport(task, drafts);
      downloadFile(`${task.name}_recommendations.md`, md, 'text/markdown');
    });
  }
  const jsonBtn = document.getElementById('showcase-export-json');
  if (jsonBtn) {
    jsonBtn.addEventListener('click', async () => {
      const drafts = await API.getDesignDrafts(taskId);
      const data = await API.getDesignJob(taskId);
      downloadFile(`${data.task.name}_recommendations.json`, JSON.stringify({ task: data.task, recommendations: drafts }, null, 2), 'application/json');
    });
  }
}

function downloadFile(filename, content, mime) {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function generateMarkdownReport(task, drafts) {
  const lines = [];
  lines.push(`# ${task.name}`);
  lines.push('');
  lines.push(`**Source:** ${task.source_job_name}`);
  lines.push(`**Created:** ${task.created_at}`);
  if (task.design_goals && task.design_goals.length) {
    lines.push(`**Goals:** ${task.design_goals.join(', ')}`);
  }
  lines.push('');
  lines.push(`## ${drafts.length} Recommendations`);
  lines.push('');
  drafts.forEach((d, i) => {
    lines.push(`### #${i + 1} · ${d.residue} (${d.chain} ${d.region})`);
    lines.push('');
    lines.push(`- **Mutation:** ${d.current_aa} → ${(d.suggested_mutations || []).join(', ')}`);
    lines.push(`- **Priority:** ${d.priority} · **Confidence:** ${d.confidence}`);
    lines.push('');
    lines.push(`**Rationale:** ${d.rationale || '-'}`);
    lines.push('');
    if (d.expected_effects && Object.keys(d.expected_effects).length) {
      lines.push('**Expected Effects:**');
      Object.entries(d.expected_effects).forEach(([k, v]) => {
        lines.push(`- ${k}: ${v}`);
      });
      lines.push('');
    }
    if (d.risks && d.risks.length) {
      lines.push('**Risks:**');
      d.risks.forEach(r => lines.push(`- ${r}`));
      lines.push('');
    }
    if (d.validation_experiments && d.validation_experiments.length) {
      lines.push('**Validation:**');
      d.validation_experiments.forEach(v => lines.push(`- ${v}`));
      lines.push('');
    }
    if (d.supporting_evidence && d.supporting_evidence.length) {
      lines.push('**Evidence:**');
      d.supporting_evidence.forEach(e => {
        if (e.type === 'md_data') {
          lines.push(`- MD: ${e.data}`);
        } else if (e.type === 'literature') {
          lines.push(`- [PMID ${e.pmid}] ${e.title}${e.relevance ? ` — ${e.relevance}` : ''}`);
        }
      });
      lines.push('');
    }
    lines.push('---');
    lines.push('');
  });
  return lines.join('\n');
}


const PROVIDER_DEFAULT_MODELS = {
  deepseek: 'deepseek-chat',
  anthropic: 'claude-opus-4-7',
  openai: 'gpt-5.5',
  gemini: 'gemini-3.1-pro',
};
const PROVIDER_DEFAULT_BASE_URLS = {
  deepseek: 'https://api.deepseek.com/v1',
  anthropic: '',
  openai: '',
  gemini: 'https://generativelanguage.googleapis.com/v1beta/openai/',
};

async function renderSettings() {
  let status;
  try {
    status = await API.getLlmSettings();
  } catch (e) {
    return `
      <div class="panel">
        <h2>Settings</h2>
        <p class="muted" style="color: var(--danger, #b91c1c);">Failed to load settings: ${escapeHtml(e.message)}</p>
      </div>
    `;
  }
  const active = status.active || {};
  const file = status.file || {};
  const env = status.env || {};

  const sourceLabel = {
    file: '~/.immunoscope/config.json',
    env_explicit: 'env var IMMUNOSCOPE_LLM_API_KEY',
    env_provider: `env var ${active.provider?.toUpperCase()}_API_KEY`,
    none: '(none configured)',
  }[active.key_source] || active.key_source;

  const keyDot = active.has_key
    ? '<span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:#10b981;margin-right:6px;"></span>'
    : '<span style="display:inline-block;width:8px;height:8px;border-radius:50%;background:#ef4444;margin-right:6px;"></span>';

  const fileProvider = file.provider || 'deepseek';
  const fileModel = file.model || '';
  const fileBaseUrl = file.base_url || '';
  const fileMaxTokens = file.max_tokens || 4096;
  const fileTemperature = typeof file.temperature === 'number' ? file.temperature : 0.7;

  return `
    <div class="panel" style="max-width: 760px;">
      <h2>LLM Configuration</h2>
      <p class="muted">
        ${keyDot}Active: <strong>${escapeHtml(active.provider || '-')}</strong>
        · model <code>${escapeHtml(active.model || '-')}</code>
        · key from <em>${escapeHtml(sourceLabel)}</em>${active.key_last4 ? ` (••••${escapeHtml(active.key_last4)})` : ''}
      </p>

      <form id="llm-settings-form" style="display:grid; gap:14px; margin-top:16px;">
        <label>
          Provider
          <select name="provider">
            <option value="deepseek"${fileProvider === 'deepseek' ? ' selected' : ''}>DeepSeek</option>
            <option value="anthropic"${fileProvider === 'anthropic' ? ' selected' : ''}>Anthropic</option>
            <option value="openai"${fileProvider === 'openai' ? ' selected' : ''}>OpenAI</option>
            <option value="gemini"${fileProvider === 'gemini' ? ' selected' : ''}>Google Gemini</option>
          </select>
        </label>

        <label>
          Model
          <input name="model" value="${escapeHtml(fileModel)}" placeholder="${escapeHtml(PROVIDER_DEFAULT_MODELS[fileProvider] || '')}">
          <span class="muted" style="font-size:11px;">Leave empty to use provider default.</span>
        </label>

        <label>
          Base URL
          <input name="base_url" value="${escapeHtml(fileBaseUrl)}" placeholder="${escapeHtml(PROVIDER_DEFAULT_BASE_URLS[fileProvider] || '')}">
          <span class="muted" style="font-size:11px;">Leave empty for the provider's official endpoint.</span>
        </label>

        <label>
          API key
          <div style="display:flex; gap:8px; align-items:center;">
            <input id="api-key-input" name="api_key" type="password" autocomplete="off" placeholder="${file.has_file_key ? '••••' + escapeHtml(active.key_last4 || '') + ' (stored)' : 'sk-... (leave empty to use env var)'}" style="flex:1;">
            <button type="button" id="api-key-clear" class="button secondary" title="Clear stored key (fall back to env)">Clear</button>
          </div>
          <span class="muted" style="font-size:11px;">Stored in <code>~/.immunoscope/config.json</code> (chmod 600). Empty value = fall back to env var.</span>
        </label>

        <div class="grid two-col">
          <label>
            Max tokens
            <input name="max_tokens" type="number" min="128" max="32768" value="${escapeHtml(String(fileMaxTokens))}">
          </label>
          <label>
            Temperature
            <input name="temperature" type="number" step="0.1" min="0" max="2" value="${escapeHtml(String(fileTemperature))}">
          </label>
        </div>

        <div style="display:flex; gap:10px; align-items:center;">
          <button type="submit" class="button">Save</button>
          <button type="button" id="llm-test-btn" class="button secondary">Test connection</button>
          <span id="llm-save-status" class="muted"></span>
        </div>

        <div id="llm-test-result" style="display:none; padding:10px 12px; border-radius:8px; font-size:13px;"></div>
      </form>

      <hr style="margin:24px 0; border:none; border-top:1px solid var(--line);">

      <h3 style="margin-top:0;">Environment fallback</h3>
      <p class="muted" style="font-size:12px;">These environment variables are used only when the corresponding field above is empty.</p>
      <ul style="font-size:13px; line-height:1.7; list-style:none; padding:0;">
        <li>${env.explicit_set ? '✓' : '✗'} <code>IMMUNOSCOPE_LLM_API_KEY</code> (overrides all providers)</li>
        <li>${env.deepseek ? '✓' : '✗'} <code>DEEPSEEK_API_KEY</code></li>
        <li>${env.anthropic ? '✓' : '✗'} <code>ANTHROPIC_API_KEY</code></li>
        <li>${env.openai ? '✓' : '✗'} <code>OPENAI_API_KEY</code></li>
        <li>${env.gemini ? '✓' : '✗'} <code>GEMINI_API_KEY</code></li>
        <li>provider override: <code>${escapeHtml(env.provider_override || '(unset)')}</code></li>
        <li>model override: <code>${escapeHtml(env.model_override || '(unset)')}</code></li>
        <li>base_url override: <code>${escapeHtml(env.base_url_override || '(unset)')}</code></li>
      </ul>

      <p class="muted" style="font-size:11px; margin-top:16px;">Config file: <code>${escapeHtml(status.config_file)}</code> (${status.config_exists ? 'exists' : 'not yet created'})</p>
    </div>
  `;
}

function mountSettings() {
  const form = document.getElementById('llm-settings-form');
  if (!form) return;
  const providerSel = form.querySelector('select[name="provider"]');
  const modelInput = form.querySelector('input[name="model"]');
  const baseUrlInput = form.querySelector('input[name="base_url"]');
  const apiKeyInput = document.getElementById('api-key-input');
  const apiKeyClear = document.getElementById('api-key-clear');
  const saveStatus = document.getElementById('llm-save-status');
  const testBtn = document.getElementById('llm-test-btn');
  const testResult = document.getElementById('llm-test-result');

  providerSel.addEventListener('change', () => {
    const p = providerSel.value;
    modelInput.placeholder = PROVIDER_DEFAULT_MODELS[p] || '';
    baseUrlInput.placeholder = PROVIDER_DEFAULT_BASE_URLS[p] || '';
  });

  apiKeyClear.addEventListener('click', () => {
    apiKeyInput.value = '';
    apiKeyInput.dataset.cleared = 'true';
    apiKeyInput.placeholder = '(cleared — will use env on save)';
  });

  form.addEventListener('submit', async event => {
    event.preventDefault();
    saveStatus.textContent = 'Saving...';
    const payload = {
      provider: providerSel.value,
      model: modelInput.value,
      base_url: baseUrlInput.value,
      max_tokens: parseInt(form.querySelector('input[name="max_tokens"]').value, 10),
      temperature: parseFloat(form.querySelector('input[name="temperature"]').value),
    };
    if (apiKeyInput.value) {
      payload.api_key = apiKeyInput.value;
    } else if (apiKeyInput.dataset.cleared === 'true') {
      payload.api_key = '';
    }
    try {
      await API.updateLlmSettings(payload);
      saveStatus.textContent = `Saved at ${new Date().toLocaleTimeString()}`;
      saveStatus.style.color = 'var(--accent-strong, #0d9488)';
      apiKeyInput.value = '';
      delete apiKeyInput.dataset.cleared;
      setTimeout(async () => {
        const newHtml = await renderSettings();
        const main = document.getElementById('page-content');
        if (main) {
          main.innerHTML = newHtml;
          mountSettings();
        }
      }, 600);
    } catch (e) {
      saveStatus.textContent = `Failed: ${e.message}`;
      saveStatus.style.color = '#b91c1c';
    }
  });

  testBtn.addEventListener('click', async () => {
    testBtn.disabled = true;
    const original = testBtn.textContent;
    testBtn.textContent = 'Testing...';
    testResult.style.display = 'block';
    testResult.style.background = '#f1f5f9';
    testResult.style.border = '1px solid var(--line)';
    testResult.style.color = 'var(--ink)';
    testResult.textContent = 'Sending a minimal probe request to the configured LLM...';
    try {
      const result = await API.testLlmConnection();
      if (result.ok) {
        testResult.style.background = '#ecfdf5';
        testResult.style.border = '1px solid #10b981';
        testResult.style.color = '#065f46';
        testResult.innerHTML = `✓ <strong>Connected</strong> — ${escapeHtml(result.provider)} / <code>${escapeHtml(result.model)}</code> in ${result.latency_ms} ms<br>Response: <code>${escapeHtml(result.response || '(empty)')}</code>`;
      } else {
        testResult.style.background = '#fef2f2';
        testResult.style.border = '1px solid #ef4444';
        testResult.style.color = '#991b1b';
        testResult.innerHTML = `✗ <strong>Failed</strong> — ${escapeHtml(result.provider)} / <code>${escapeHtml(result.model)}</code><br>${escapeHtml(result.error || 'Unknown error')}`;
      }
    } catch (e) {
      testResult.style.background = '#fef2f2';
      testResult.style.border = '1px solid #ef4444';
      testResult.style.color = '#991b1b';
      testResult.innerHTML = `✗ <strong>Request failed</strong> — ${escapeHtml(e.message)}`;
    } finally {
      testBtn.disabled = false;
      testBtn.textContent = original;
    }
  });
}

function renderDocs() {
  return `
    <div class="panel">
      <h2>Documentation</h2>
      <p class="muted">Learn how to use ImmunoScope.</p>
      <p>This feature will be integrated soon.</p>
    </div>
  `;
}

// Agent WebSocket Client
class AgentClient {
  constructor() {
    this.ws = null;
    this.sessionId = null;
    this.isConnected = false;
    this.reconnectAttempts = 0;
    this.maxReconnectAttempts = 5;
    this.currentMessageDiv = null;
    this.init();
  }

  init() {
    this.connectWebSocket();
    this.setupEventListeners();
  }

  connectWebSocket() {
    const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
    const wsUrl = `${protocol}//${window.location.host}/api/agent/ws`;

    this.ws = new WebSocket(wsUrl);

    this.ws.onopen = () => {
      console.log('WebSocket connected');
      this.isConnected = true;
      this.reconnectAttempts = 0;
      this.updateStatus('Connected', true);
    };

    this.ws.onmessage = (event) => {
      try {
        const message = JSON.parse(event.data);
        this.handleMessage(message);
      } catch (error) {
        console.error('Failed to parse message:', error);
      }
    };

    this.ws.onerror = (error) => {
      console.error('WebSocket error:', error);
      this.showError('Connection error occurred');
    };

    this.ws.onclose = () => {
      console.log('WebSocket disconnected');
      this.isConnected = false;
      this.updateStatus('Disconnected', false);

      if (this.reconnectAttempts < this.maxReconnectAttempts) {
        this.reconnectAttempts++;
        const delay = Math.min(1000 * Math.pow(2, this.reconnectAttempts), 30000);
        console.log(`Reconnecting in ${delay}ms (attempt ${this.reconnectAttempts})`);
        setTimeout(() => this.connectWebSocket(), delay);
      } else {
        this.showError('Failed to reconnect. Please refresh the page.');
      }
    };
  }

  handleMessage(message) {
    switch (message.type) {
      case 'session_info':
        this.sessionId = message.session_id;
        console.log('Session ID:', this.sessionId);
        break;

      case 'text_delta':
        this.appendAgentText(message.text);
        break;

      case 'tool_use_start':
        this.showToolCall(message.name, message.id, 'running');
        break;

      case 'tool_use_end':
        const status = message.ok ? 'success' : 'error';
        this.updateToolCall(message.id, status, message.name);
        break;

      case 'turn_complete':
        this.finishAgentMessage();
        this.hideTypingIndicator();
        this.enableInput();
        break;

      case 'error':
        this.showError(message.message || 'An error occurred');
        this.hideTypingIndicator();
        this.enableInput();
        break;

      case 'aborted':
        this.appendAgentText('\n\n[Operation aborted by user]');
        this.finishAgentMessage();
        this.hideTypingIndicator();
        this.enableInput();
        break;

      case 'pong':
        break;

      default:
        console.warn('Unknown message type:', message.type);
    }
  }

  sendMessage(content) {
    if (!this.isConnected) {
      this.showError('Not connected to server');
      return;
    }

    if (!content.trim()) {
      return;
    }

    this.addUserMessage(content);

    this.ws.send(JSON.stringify({
      type: 'message',
      content: content
    }));

    this.showTypingIndicator();
    this.disableInput();

    const input = document.getElementById('messageInput');
    if (input) {
      input.value = '';
      this.adjustTextareaHeight();
    }
  }

  addUserMessage(content) {
    const chatContainer = document.getElementById('chatContainer');

    const welcome = chatContainer.querySelector('.welcome-message');
    if (welcome) {
      welcome.remove();
    }

    const messageDiv = document.createElement('div');
    messageDiv.className = 'message user';
    messageDiv.innerHTML = `
      <div class="message-content">${this.escapeHtml(content)}</div>
      <div class="message-avatar">U</div>
    `;

    chatContainer.appendChild(messageDiv);
    this.scrollToBottom();
  }

  startAgentMessage() {
    if (this.currentMessageDiv) {
      return;
    }

    const chatContainer = document.getElementById('chatContainer');

    this.currentMessageDiv = document.createElement('div');
    this.currentMessageDiv.className = 'message agent';
    this.currentMessageDiv.innerHTML = `
      <div class="message-avatar">A</div>
      <div class="message-content"></div>
    `;

    chatContainer.appendChild(this.currentMessageDiv);
    this.scrollToBottom();
  }

  appendAgentText(text) {
    if (!this.currentMessageDiv) {
      this.startAgentMessage();
    }

    const contentDiv = this.currentMessageDiv.querySelector('.message-content');
    contentDiv.textContent += text;
    this.scrollToBottom();
  }

  finishAgentMessage() {
    // End of turn: close out the "Working" fold (option A) and clear the bubble.
    this._finalizeProcess();
    this.currentMessageDiv = null;
  }

  // ---- Agent "process" fold (option A), mirroring the ImmunoScope design chat ----
  // Tool calls + intermediate narration collapse into one muted strip so the
  // final reply bubble reads clean.
  _getProcessBlock() {
    const chatContainer = document.getElementById('chatContainer');
    if (this.processBlock && chatContainer.contains(this.processBlock)) return this.processBlock;
    if (typeof _maEnsureSpinStyle === 'function') _maEnsureSpinStyle();
    const welcome = chatContainer.querySelector('.welcome-message');
    if (welcome) welcome.remove();
    const wrap = document.createElement('div');
    wrap.className = 'agent-process';
    wrap.style.cssText = 'padding: 2px 16px 2px 60px;';
    wrap.innerHTML = `
      <details style="border:1px solid var(--line); border-radius:8px; background:#fafafa;">
        <summary style="cursor:pointer; list-style:none; display:flex; align-items:center; gap:7px; padding:6px 10px; font-size:12px; color:var(--muted); user-select:none;">
          <span class="ap-spin" style="display:inline-block; width:12px; height:12px; border:2px solid var(--muted); border-top-color:transparent; border-radius:50%; animation:maspin 0.8s linear infinite;"></span>
          <span>Working · <span class="ap-count">0</span> steps</span>
          <span style="margin-left:auto; font-size:10px; opacity:0.7;">expand</span>
        </summary>
        <div class="ap-body" style="padding:6px 12px 8px 28px; font-family: ui-monospace, monospace; font-size:11px; color:#6f746d; line-height:1.7; border-top:1px solid var(--line); white-space:pre-wrap; word-break:break-word;"></div>
      </details>`;
    chatContainer.appendChild(wrap);
    this.processBlock = wrap;
    this.processStepCount = 0;
    return wrap;
  }

  _addProcessStep(label, toolId) {
    const block = this._getProcessBlock();
    const line = document.createElement('div');
    line.textContent = label;
    if (toolId) line.dataset.toolId = toolId;
    block.querySelector('.ap-body').appendChild(line);
    this.processStepCount = (this.processStepCount || 0) + 1;
    const cnt = block.querySelector('.ap-count');
    if (cnt) cnt.textContent = String(this.processStepCount);
    this.scrollToBottom();
    return line;
  }

  _foldNarration() {
    if (!this.currentMessageDiv) return;
    const contentDiv = this.currentMessageDiv.querySelector('.message-content');
    const txt = (contentDiv ? contentDiv.textContent : '').trim();
    if (txt) this._addProcessStep('💬 ' + txt);
    this.currentMessageDiv.remove();
    this.currentMessageDiv = null;
  }

  _finalizeProcess() {
    if (!this.processBlock) return;
    const spin = this.processBlock.querySelector('.ap-spin');
    if (spin) {
      spin.style.animation = 'none';
      spin.style.border = 'none';
      spin.style.width = 'auto';
      spin.style.height = 'auto';
      spin.style.color = '#255b4b';
      spin.style.fontWeight = '800';
      spin.textContent = '✓';
    }
    this.processBlock = null;
    this.processStepCount = 0;
  }

  showToolCall(toolName, toolId, status = 'running') {
    // Fold any in-progress narration + this tool call into the muted strip.
    this._foldNarration();
    this._addProcessStep('⚙ ' + toolName, toolId);
  }

  updateToolCall(toolId, status, toolName) {
    if (status !== 'error' || !this.processBlock) return;
    const body = this.processBlock.querySelector('.ap-body');
    if (!body) return;
    body.querySelectorAll('div').forEach(line => {
      if (line.dataset && line.dataset.toolId === toolId && !/✗$/.test(line.textContent)) {
        line.textContent += '  ✗';
      }
    });
  }

  showTypingIndicator() {
    const chatContainer = document.getElementById('chatContainer');

    let indicator = document.getElementById('typingIndicator');
    if (!indicator) {
      indicator = document.createElement('div');
      indicator.id = 'typingIndicator';
      indicator.className = 'message agent';
      indicator.innerHTML = `
        <div class="message-avatar">A</div>
        <div class="typing-indicator active">
          <div class="typing-dots">
            <span></span>
            <span></span>
            <span></span>
          </div>
        </div>
      `;
      chatContainer.appendChild(indicator);
    }

    indicator.style.display = 'flex';
    this.scrollToBottom();
  }

  hideTypingIndicator() {
    const indicator = document.getElementById('typingIndicator');
    if (indicator) {
      indicator.remove();
    }
  }

  updateStatus(text, connected) {
    const statusText = document.getElementById('statusText');
    const statusDot = document.getElementById('statusDot');

    if (statusText) statusText.textContent = text;

    if (statusDot) {
      if (connected) {
        statusDot.classList.remove('disconnected');
      } else {
        statusDot.classList.add('disconnected');
      }
    }
  }

  showError(message) {
    const errorBanner = document.getElementById('errorBanner');
    if (errorBanner) {
      errorBanner.textContent = message;
      errorBanner.classList.add('active');
      setTimeout(() => this.hideError(), 5000);
    }
  }

  hideError() {
    const errorBanner = document.getElementById('errorBanner');
    if (errorBanner) {
      errorBanner.classList.remove('active');
    }
  }

  enableInput() {
    const input = document.getElementById('messageInput');
    const button = document.getElementById('sendButton');
    if (input) input.disabled = false;
    if (button) button.disabled = false;
  }

  disableInput() {
    const input = document.getElementById('messageInput');
    const button = document.getElementById('sendButton');
    if (input) input.disabled = true;
    if (button) button.disabled = true;
  }

  scrollToBottom() {
    const chatContainer = document.getElementById('chatContainer');
    if (chatContainer) {
      chatContainer.scrollTop = chatContainer.scrollHeight;
    }
  }

  adjustTextareaHeight() {
    const textarea = document.getElementById('messageInput');
    if (textarea) {
      textarea.style.height = 'auto';
      textarea.style.height = Math.min(textarea.scrollHeight, 150) + 'px';
    }
  }

  escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
  }

  setupEventListeners() {
    const input = document.getElementById('messageInput');
    const button = document.getElementById('sendButton');

    if (input) {
      input.addEventListener('input', () => this.adjustTextareaHeight());

      input.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' && !e.shiftKey) {
          e.preventDefault();
          this.sendMessage(input.value);
        }
      });
    }

    if (button) {
      button.addEventListener('click', () => {
        this.sendMessage(input.value);
      });
    }

    setInterval(() => {
      if (this.isConnected) {
        this.ws.send(JSON.stringify({ type: 'ping' }));
      }
    }, 30000);
  }

  saveHistory() {
    const chatContainer = document.getElementById('chatContainer');
    if (!chatContainer) return null;

    // Save the entire chat HTML
    return {
      html: chatContainer.innerHTML,
      sessionId: this.sessionId
    };
  }

  restoreHistory(history) {
    if (!history) return;

    const chatContainer = document.getElementById('chatContainer');
    if (!chatContainer) return;

    // Restore chat HTML
    chatContainer.innerHTML = history.html;
    this.sessionId = history.sessionId;

    // Scroll to bottom
    this.scrollToBottom();
  }

  cleanup() {
    if (this.ws) {
      this.ws.close();
      this.ws = null;
    }
    this.isConnected = false;
    // Don't clear sessionId - keep it for history restoration
    this.currentMessageDiv = null;
  }
}

window.App = App;
document.addEventListener('DOMContentLoaded', () => App.start());
