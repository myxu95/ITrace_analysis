class ReportsManager {
    constructor() {
        this.reports = [];
        this.filteredReports = [];
        this.currentFilter = 'all';
        this.searchQuery = '';

        this.init();
    }

    init() {
        this.setupEventListeners();
        this.loadReports();
    }

    setupEventListeners() {
        document.getElementById('refreshBtn').addEventListener('click', () => {
            this.loadReports();
        });

        document.getElementById('searchInput').addEventListener('input', (e) => {
            this.searchQuery = e.target.value.toLowerCase();
            this.filterReports();
        });

        document.querySelectorAll('.filter-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
                e.target.classList.add('active');
                this.currentFilter = e.target.dataset.filter;
                this.filterReports();
            });
        });

        document.getElementById('closeViewer').addEventListener('click', () => {
            this.closeViewer();
        });

        document.getElementById('viewerOverlay').addEventListener('click', (e) => {
            if (e.target.id === 'viewerOverlay') {
                this.closeViewer();
            }
        });

        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                this.closeViewer();
            }
        });
    }

    async loadReports() {
        const container = document.getElementById('reportsContainer');
        container.innerHTML = '<div class="loading">Loading reports...</div>';

        try {
            const response = await fetch('/api/reports/list');
            if (!response.ok) {
                throw new Error(`HTTP ${response.status}`);
            }

            const data = await response.json();
            this.reports = data.reports || [];
            this.filterReports();
        } catch (error) {
            console.error('Failed to load reports:', error);
            container.innerHTML = `
                <div class="empty-state">
                    <h2>⚠️ Failed to Load Reports</h2>
                    <p>${error.message}</p>
                </div>
            `;
        }
    }

    filterReports() {
        this.filteredReports = this.reports.filter(report => {
            const matchesFilter = this.currentFilter === 'all' || report.type === this.currentFilter;
            const matchesSearch = !this.searchQuery ||
                report.name.toLowerCase().includes(this.searchQuery) ||
                report.path.toLowerCase().includes(this.searchQuery);
            return matchesFilter && matchesSearch;
        });

        this.renderReports();
    }

    renderReports() {
        const container = document.getElementById('reportsContainer');

        if (this.filteredReports.length === 0) {
            container.innerHTML = `
                <div class="empty-state">
                    <h2>📭 No Reports Found</h2>
                    <p>No reports match your current filters.</p>
                </div>
            `;
            return;
        }

        const grid = document.createElement('div');
        grid.className = 'reports-grid';

        this.filteredReports.forEach(report => {
            const card = this.createReportCard(report);
            grid.appendChild(card);
        });

        container.innerHTML = '';
        container.appendChild(grid);
    }

    createReportCard(report) {
        const card = document.createElement('div');
        card.className = 'report-card';
        card.onclick = () => this.openReport(report);

        const typeClass = report.type === 'comparison' ? 'comparison' : 'analysis';
        const typeLabel = report.type === 'comparison' ? 'Comparison' : 'Analysis';

        card.innerHTML = `
            <span class="type-badge ${typeClass}">${typeLabel}</span>
            <h3>${this.escapeHtml(report.name)}</h3>
            <div class="meta">
                <div class="meta-item">
                    <span>📅</span>
                    <span>${this.formatDate(report.modified)}</span>
                </div>
                <div class="meta-item">
                    <span>📦</span>
                    <span>${this.formatSize(report.size)}</span>
                </div>
            </div>
            <div class="path" title="${this.escapeHtml(report.path)}">${this.escapeHtml(report.relative_path)}</div>
        `;

        return card;
    }

    openReport(report) {
        const overlay = document.getElementById('viewerOverlay');
        const iframe = document.getElementById('viewerIframe');
        const title = document.getElementById('viewerTitle');

        title.textContent = report.name;
        iframe.src = `/api/reports/view?path=${encodeURIComponent(report.path)}`;
        overlay.classList.add('active');
    }

    closeViewer() {
        const overlay = document.getElementById('viewerOverlay');
        const iframe = document.getElementById('viewerIframe');

        overlay.classList.remove('active');
        iframe.src = '';
    }

    formatDate(timestamp) {
        const date = new Date(timestamp * 1000);
        const now = new Date();
        const diffMs = now - date;
        const diffDays = Math.floor(diffMs / (1000 * 60 * 60 * 24));

        if (diffDays === 0) {
            return 'Today';
        } else if (diffDays === 1) {
            return 'Yesterday';
        } else if (diffDays < 7) {
            return `${diffDays} days ago`;
        } else {
            return date.toLocaleDateString();
        }
    }

    formatSize(bytes) {
        if (bytes < 1024) return bytes + ' B';
        if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
        return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
    }

    escapeHtml(text) {
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }
}

document.addEventListener('DOMContentLoaded', () => {
    new ReportsManager();
});
