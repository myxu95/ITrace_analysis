import React from 'react';
import {
  LayoutDashboard,
  FolderKanban,
  Activity,
  FlaskConical,
  GitCompareArrows,
  FileText,
  Waves,
  Sparkles,
  BrainCircuit,
  Database,
  Settings,
  BookOpen,
  ChevronRight,
} from 'lucide-react';

const Sidebar = () => {
  const [activeItem, setActiveItem] = React.useState('dashboard');

  const navigationSections = [
    {
      title: 'WORKSPACE',
      items: [
        { id: 'dashboard', label: 'Dashboard', icon: LayoutDashboard, href: '#/' },
        { id: 'projects', label: 'Projects', icon: FolderKanban, href: '#/projects' },
        { id: 'jobs', label: 'Jobs', icon: Activity, href: '#/jobs', badge: '3' }
      ]
    },
    {
      title: 'ANALYSIS',
      items: [
        { id: 'new', label: 'New Analysis', icon: FlaskConical, href: '#/new' },
        { id: 'compare', label: 'Compare', icon: GitCompareArrows, href: '#/compare' },
        { id: 'reports', label: 'Reports', icon: FileText, href: '#/reports' },
        { id: 'trajectories', label: 'Trajectories', icon: Waves, href: '#/trajectories' }
      ]
    },
    {
      title: 'AI',
      items: [
        { id: 'agent', label: 'Ask ImmunoScope', icon: Sparkles, href: '#/agent', highlight: true },
        { id: 'copilot', label: 'Design Copilot', icon: BrainCircuit, href: '#/copilot' },
        { id: 'knowledge', label: 'Knowledge Base', icon: Database, href: '#/knowledge' }
      ]
    },
    {
      title: 'SYSTEM',
      items: [
        { id: 'settings', label: 'Settings', icon: Settings, href: '#/settings' },
        { id: 'docs', label: 'Documentation', icon: BookOpen, href: '#/docs' }
      ]
    }
  ];

  const currentWorkspace = {
    name: 'A0201_TAX_JM22',
    trajectories: 3,
    duration: '200 ns × 3 replicas',
    status: 'Analysis completed',
    statusColor: 'text-teal-600'
  };

  return (
    <div className="flex flex-col h-screen w-[260px] bg-white border-r border-slate-200">
      {/* Header */}
      <div className="flex items-center gap-3 px-5 py-6 border-b border-slate-100">
        <div className="w-10 h-10 rounded-full bg-white flex items-center justify-center shadow-sm ring-1 ring-slate-200 overflow-hidden">
          <img
            src="/static/immunoscope-logo-icon.png"
            alt="ImmunoScope logo"
            className="w-full h-full object-cover"
          />
        </div>
        <div className="flex flex-col">
          <h1 className="text-base font-semibold text-slate-900 tracking-tight">
            ImmunoScope
          </h1>
          <p className="text-xs text-slate-500 font-medium">
            MD Analysis Platform
          </p>
        </div>
      </div>

      {/* Navigation */}
      <nav className="flex-1 overflow-y-auto px-3 py-4 space-y-6">
        {navigationSections.map((section) => (
          <div key={section.title}>
            <h2 className="px-3 mb-2 text-[10px] font-bold text-slate-400 uppercase tracking-wider">
              {section.title}
            </h2>
            <div className="space-y-0.5">
              {section.items.map((item) => {
                const Icon = item.icon;
                const isActive = activeItem === item.id;

                return (
                  <a
                    key={item.id}
                    href={item.href}
                    onClick={() => setActiveItem(item.id)}
                    className={`
                      group relative flex items-center gap-3 px-3 py-2 rounded-xl
                      transition-all duration-200 ease-out
                      ${isActive
                        ? 'bg-teal-50 text-teal-700'
                        : 'text-slate-600 hover:bg-slate-50 hover:text-slate-900'
                      }
                    `}
                  >
                    {/* Active indicator */}
                    {isActive && (
                      <div className="absolute left-0 top-1/2 -translate-y-1/2 w-0.5 h-5 bg-teal-600 rounded-r-full" />
                    )}

                    <Icon
                      className={`
                        w-[18px] h-[18px] flex-shrink-0
                        transition-colors duration-200
                        ${isActive ? 'text-teal-600' : 'text-slate-400 group-hover:text-slate-600'}
                      `}
                      strokeWidth={2}
                    />

                    <span className={`
                      text-[13px] font-medium flex-1
                      ${isActive ? 'font-semibold' : ''}
                    `}>
                      {item.label}
                    </span>

                    {/* Highlight badge for AI features */}
                    {item.highlight && !isActive && (
                      <div className="w-1.5 h-1.5 rounded-full bg-teal-500 animate-pulse" />
                    )}

                    {/* Badge */}
                    {item.badge && (
                      <span className="px-1.5 py-0.5 text-[10px] font-bold bg-slate-100 text-slate-600 rounded-md">
                        {item.badge}
                      </span>
                    )}

                    {/* Hover arrow */}
                    <ChevronRight
                      className={`
                        w-3.5 h-3.5 text-slate-300 opacity-0 group-hover:opacity-100
                        transition-opacity duration-200
                        ${isActive ? 'hidden' : ''}
                      `}
                    />
                  </a>
                );
              })}
            </div>
          </div>
        ))}
      </nav>

      {/* Current Workspace Card */}
      <div className="px-3 pb-4">
        <div className="bg-gradient-to-br from-slate-50 to-slate-100/50 rounded-xl p-4 border border-slate-200/60 shadow-sm">
          <div className="flex items-start justify-between mb-3">
            <div>
              <p className="text-[10px] font-bold text-slate-400 uppercase tracking-wider mb-1">
                Current Workspace
              </p>
              <h3 className="text-sm font-semibold text-slate-900 font-mono">
                {currentWorkspace.name}
              </h3>
            </div>
            <div className="w-2 h-2 rounded-full bg-teal-500 animate-pulse" />
          </div>

          <div className="space-y-2 mb-3">
            <div className="flex items-center gap-2 text-xs text-slate-600">
              <Waves className="w-3.5 h-3.5 text-slate-400" strokeWidth={2} />
              <span className="font-medium">{currentWorkspace.trajectories} trajectories</span>
            </div>
            <div className="text-xs text-slate-500 font-mono pl-5">
              {currentWorkspace.duration}
            </div>
          </div>

          <div className="pt-3 border-t border-slate-200/60">
            <div className="flex items-center justify-between">
              <span className="text-[11px] font-medium text-slate-500">Status</span>
              <span className={`text-[11px] font-semibold ${currentWorkspace.statusColor}`}>
                {currentWorkspace.status}
              </span>
            </div>
          </div>
        </div>
      </div>

      {/* Footer */}
      <div className="px-5 py-3 border-t border-slate-100">
        <div className="flex items-center gap-2">
          <div className="w-7 h-7 rounded-full bg-gradient-to-br from-slate-200 to-slate-300 flex items-center justify-center">
            <span className="text-[11px] font-bold text-slate-600">XM</span>
          </div>
          <div className="flex-1 min-w-0">
            <p className="text-xs font-medium text-slate-700 truncate">xumy</p>
            <p className="text-[10px] text-slate-400">Local session</p>
          </div>
        </div>
      </div>
    </div>
  );
};

export default Sidebar;
