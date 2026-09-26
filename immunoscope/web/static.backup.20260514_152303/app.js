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
  getJob(id) {
    return this.request('GET', `/jobs/${id}`);
  },
  createJob(formData) {
    return this.request('POST', '/jobs', formData, true);
  },
  createCompareJob(payload) {
    return this.request('POST', '/compare/jobs', payload);
  },
  startJob(id) {
    return this.request('POST', `/jobs/${id}/start`);
  },
  deleteJob(id) {
    return this.request('DELETE', `/jobs/${id}`);
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
  standard: 'quality,identity,contact,rmsf,bsa,rrcs,cluster,landscape,report',
  contact: 'identity,contact,rrcs,report',
  fast: 'quality,identity,contact,report',
  custom: 'quality,identity,contact,cluster,report',
};

const QUICK_QUESTIONS = [
  'Summarize this trajectory',
  'Diagnose stability',
  'Summarize interface',
  'Identify hotspots',
  'Suggest mutation sites',
];

const App = {
  async start() {
    window.addEventListener('hashchange', () => this.navigate());
    await this.navigate();
    window.setInterval(() => {
      const hash = location.hash.slice(1) || '/';
      if (hash === '/jobs' || hash.startsWith('/jobs/')) this.navigate(false);
    }, 6000);
  },
  async navigate(showLoading = true) {
    const hash = (location.hash.slice(1) || '/').replace(/\/$/, '') || '/';
    const content = document.getElementById('page-content');
    if (showLoading) content.innerHTML = '<div class="panel">Loading...</div>';
    document.querySelectorAll('.nav a').forEach(link => {
      link.classList.toggle('active', link.getAttribute('href') === `#${hash}`);
    });
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
    if (hash === '/compare') {
      document.getElementById('page-title').textContent = 'Compare';
      content.innerHTML = await renderCompare();
      mountCompare();
      return;
    }
    if (hash.startsWith('/jobs/')) {
      const id = hash.split('/')[2];
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
  return `
    <form id="compare-form" class="panel analysis-form">
      <div class="form-header">
        <div>
          <h2>Create Compare Job</h2>
          <p class="muted">Compare two existing ImmunoScope job/report roots. Use this for Standard MD vs REST2, replicate checks, WT vs mutant, or binder comparisons.</p>
        </div>
        <label class="switch-row">
          <input name="auto_start" type="checkbox" checked value="true">
          Start immediately
        </label>
      </div>
      <div class="grid two-col">
        <label>Case A root<input name="case_a" placeholder="/path/to/output/job_standard"></label>
        <label>Case B root<input name="case_b" placeholder="/path/to/output/job_rest2"></label>
      </div>
      <div class="grid two-col">
        <label>Label A<input name="label_a" value="Standard MD"></label>
        <label>Label B<input name="label_b" value="Enhanced sampling"></label>
      </div>
      <div class="grid three-col">
        <label>
          Mode
          <select name="comparison_mode">
            <option value="sampling">Sampling</option>
            <option value="replicate">Replicate</option>
            <option value="mutation">Mutation</option>
            <option value="generic">Generic</option>
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
        <label>Alignment selection<input name="alignment_selection" value="phla_core_ca"></label>
      </div>
      <label>
        Residue mapping policy
        <select name="residue_mapping">
          <option value="auto">Auto</option>
          <option value="identity">Identity</option>
          <option value="region-only">Region only</option>
        </select>
      </label>
      <label>
        Context
        <textarea name="comparison_context" placeholder="Optional comparison context shown in the report"></textarea>
      </label>
      <button type="submit">Create compare job</button>
      <p id="compare-status" class="muted"></p>
    </form>
  `;
}

function mountCompare() {
  const form = document.getElementById('compare-form');
  const status = document.getElementById('compare-status');
  form.addEventListener('submit', async event => {
    event.preventDefault();
    status.textContent = 'Creating compare job...';
    const data = new FormData(form);
    const payload = Object.fromEntries(data.entries());
    payload.auto_start = form.querySelector('input[name="auto_start"]').checked;
    try {
      const job = await API.createCompareJob(payload);
      status.textContent = `Created compare job ${job.id}.`;
      location.hash = `#/jobs/${job.id}`;
    } catch (error) {
      status.textContent = `Failed: ${error.message}`;
    }
  });
}

async function renderNewAnalysis() {
  return `
    <form id="new-job-form" class="panel analysis-form">
      <div class="form-header">
        <div>
          <h2>Create Analysis Job</h2>
          <p class="muted">Upload raw or prepared trajectory inputs. ImmunoScope will validate the inputs and execute the selected modules through the CLI workflow.</p>
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
            <option value="raw">Raw trajectory: structure + topology + trajectory</option>
            <option value="prepared">Prepared trajectory: prepared structure + topology + processed trajectory</option>
          </select>
        </label>
        <label>
          Module preset
          <select id="module-preset">
            <option value="standard">Standard full report</option>
            <option value="contact">Contact/RRCS report</option>
            <option value="fast">Fast intake report</option>
            <option value="custom">Custom modules</option>
          </select>
        </label>
      </div>

      <section class="input-block" data-mode="raw">
        <h3>Raw Inputs</h3>
        <div class="grid three-col">
          <label>Structure PDB/GRO<input name="structure" type="file" accept=".pdb,.gro"></label>
          <label>Topology TPR<input name="topology" type="file" accept=".tpr"></label>
          <label>Trajectory XTC/TRR<input name="trajectory" type="file" accept=".xtc,.trr"></label>
        </div>
      </section>

      <section class="input-block hidden" data-mode="prepared">
        <h3>Prepared Inputs</h3>
        <div class="grid three-col">
          <label>Prepared structure<input name="prepared_structure" type="file" accept=".pdb,.gro"></label>
          <label>Topology TPR<input name="topology_prepared" type="file" accept=".tpr"></label>
          <label>Processed trajectory<input name="processed_trajectory" type="file" accept=".xtc,.trr"></label>
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

      <label>
        Modules
        <input name="modules" id="modules-input" value="${MODULE_PRESETS.standard}">
      </label>
      <label>
        Notes
        <textarea name="notes" placeholder="Optional context for this job"></textarea>
      </label>
      <button type="submit">Create job</button>
      <p id="new-job-status" class="muted"></p>
    </form>
  `;
}

function mountNewAnalysis() {
  const form = document.getElementById('new-job-form');
  const status = document.getElementById('new-job-status');
  const mode = document.getElementById('input-mode');
  const preset = document.getElementById('module-preset');
  const modules = document.getElementById('modules-input');

  function updateMode() {
    document.querySelectorAll('.input-block').forEach(block => {
      block.classList.toggle('hidden', block.dataset.mode !== mode.value);
    });
  }
  mode.addEventListener('change', updateMode);
  preset.addEventListener('change', () => {
    if (preset.value !== 'custom') modules.value = MODULE_PRESETS[preset.value];
  });
  modules.addEventListener('input', () => { preset.value = 'custom'; });

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
        <h3>ImmunoScope Assistant</h3>
        <div id="assistant-skill-strip" class="assistant-skill-strip">
          <span class="muted">Loading available skills...</span>
        </div>
        <div class="quick-actions">
          ${QUICK_QUESTIONS.map(q => `<button class="button secondary quick-question" data-question="${escapeHtml(q)}">${escapeHtml(q)}</button>`).join('')}
        </div>
        <div id="assistant-messages" class="assistant-messages">
          <div class="assistant-message assistant">Ask about this job after or during analysis. Answers use available structured results and fall back gracefully when evidence is incomplete.</div>
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

document.addEventListener('DOMContentLoaded', () => App.start());
